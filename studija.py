"""DFT/DFPT monosloja TiSe2. Primjer: python studija.py run --temperatures 100 --part m"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager, nullcontext
from datetime import datetime, timedelta
import json
import os
from pathlib import Path
import re
import shutil
from threading import Event, Thread
from time import monotonic

import numpy as np
from filelock import FileLock, Timeout
from scipy.constants import atomic_mass, elementary_charge

import qe
import modovi as mod


STATUS_INTERVAL_SECONDS = 20 * 60


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def tdir(c, t):
    return qe.run_root(c) / f"T{t:04d}K"


def settings(path):
    c = read_json(path)
    ts = c["temperatures_K"]
    if not ts or list(sorted(set(ts))) != ts or min(ts) <= 0:
        raise ValueError("Temperature moraju biti pozitivne, različite i poredane.")
    k, q = c["k_mesh"], c["q_mesh"]
    if k < 2 or q < 2 or k % 2 or q % 2 or k % q:
        raise ValueError("Izaberi parne k i q mreže uz k djeljiv s q. M je tada na q-mreži.")
    a = np.array(c["amplitudes_A"])
    if len(a) < 3 or not np.all(np.diff(a) > 0) or 0 not in a or not np.allclose(a, -a[::-1]):
        raise ValueError("Za središnju razliku trebaju uređeni simetrični pomaci s nulom.")
    return c


def release_scratch(directory, c):
    # Uklanjaju se samo obnovljive QE valne funkcije nakon spremanja trajnih izlaza.
    p = directory / "tmp"
    if (
        p.is_dir()
        and not p.is_symlink()
        and p.resolve().parent == directory.resolve()
        and qe.run_root(c).resolve() in p.resolve().parents
    ):
        shutil.rmtree(p)


def relax(c, t, g, d):
    """BFGS za ionske koordinate uz fiksnu ćeliju. Spremi konačne položaje u Å."""
    final = d / "geometrija.json"
    if final.exists() and (d / "relax.ok.json").exists():
        if not (d / "energija.json").exists():
            qe.save_json(d / "energija.json", {"F_eV": qe.energy(d / "relax.out")})
            qe.xyz(d / "pocetna.extxyz", g)
            release_scratch(d, c)
        return read_json(final)
    qe.xyz(d / "pocetna.extxyz", g)
    log = qe.step(
        d,
        "relax",
        "pw",
        qe.pw_input(g, c, t, d, relax=True),
        ["tmp/tise2.save/data-file-schema.xml"],
        kind="relax",
    )
    g2 = qe.final_geometry(d, g)
    qe.xyz(d / "geometrija.extxyz", g2)
    qe.save_json(d / "energija.json", {"F_eV": qe.energy(log)})
    qe.save_json(final, g2)
    release_scratch(d, c)
    return g2


def scf(c, t, g, d):
    return qe.step(
        d, "scf", "pw", qe.pw_input(g, c, t, d),
        ["tmp/tise2.save/data-file-schema.xml"], kind="phonon_scf",
    )


def primitive_part(c, t, g, part):
    """Izračunaj najniži M izravno ili devet grana iz Fourierovih konstanti sila."""
    d = tdir(c, t) / "fononi_1x1"
    if part == "m" and (d / "M_direktno.json").exists():
        return
    if part == "mesh" and (d / "minimum.json").exists():
        return
    scf(c, t, g, d)
    if part == "m":
        log = qe.step(d, "M", "ph", qe.ph_input(c, "pointM"), ["pointM.dyn"])
        # Frekvencije nakon konačne dijagonalizacije, s predznakom prema konvenciji QE.
        pat = r"freq\s*\(\s*\d+\s*\)\s*=\s*(" + qe.REAL + r")\s*\[THz\]"
        vals = [float(s) for s in re.findall(pat, log.read_text(errors="replace"), re.I)]
        if len(vals) < 9:
            raise ValueError(f"Nedostaje devet frekvencija u {log}.")
        qe.save_json(
            d / "M_direktno.json",
            {"q_frac": [0.5, 0, 0], "nu_THz": vals[-9:], "nu_min_THz": min(vals[-9:])},
        )
        return
    qe.step(d, "mesh", "ph", qe.ph_input(c, "mesh"), ["mesh.dyn0"])
    qe.step(d, "q2r", "q2r", qe.Q2R, ["mesh.fc"])
    qe.step(d, "matdyn", "matdyn", qe.matdyn_input(c), ["path.freq"])
    qe.save_json(d / "minimum.json", mod.dispersion(d, c, g))
    release_scratch(d, c)


def gamma(c, t, g, d, reference):
    """Novi elektronski SCF i DFPT u Γ za točno onu strukturu koja se ispituje."""
    if (d / "gamma.npz").exists() and (d / "zavrseno.json").exists():
        with np.load(d / "gamma.npz") as a:
            return a["D_signed_THz2"]
    qe.xyz(d / "struktura.extxyz", g)
    log = scf(c, t, g, d)
    qe.step(d, "gamma", "ph", qe.ph_input(c, "gamma"), ["gamma.dyn"])
    qe.step(d, "dynmat", "dynmat", qe.DYNMAT, ["gamma.eig"])
    mat = mod.super_gamma(d, g, c["instability_threshold_THz"], reference)
    qe.save_json(d / "zavrseno.json", {"F_eV": qe.energy(log)})
    release_scratch(d, c)
    return mat


def reference(c, t, g):
    d = tdir(c, t) / "gamma_2x2_ref"
    return gamma(c, t, qe.supercell(g), d, True)


def relaxed(c, t, g):
    g4 = qe.supercell(g)
    g2 = relax(c, t, qe.distorted(g4, c), tdir(c, t) / "relax_2x2")
    gamma(c, t, g2, tdir(c, t) / "gamma_2x2_rel", False)


def directions(c, t, g):
    """Odredi predznake na prvoj T, zatim prati svaku granu u njezinu sektoru."""
    ts = c["temperatures_K"]
    folder = tdir(c, t) / "gamma_2x2_ref"
    previous = None
    if t != ts[0]:
        prev = ts[ts.index(t) - 1]
        prev_file = tdir(c, prev) / "gamma_2x2_ref" / "smjerovi_M.json"
        if not prev_file.exists():
            with lock(c, prev):
                if not prev_file.exists():
                    gprev = relax(c, prev, qe.primitive(c), tdir(c, prev) / "relax_1x1")
                    directions(c, prev, gprev)
        previous = read_json(prev_file)
    g4 = qe.supercell(g)
    info = mod.frozen_directions(reference(c, t, g), g4, c["instability_threshold_THz"], previous)
    for name, item in info["smjerovi"].items():
        vec = np.asarray(item["e_massweighted"]).reshape(-1, 3)
        unit_q = vec / np.sqrt(mod.masses(g4))[:, None]
        mod.csv_out(
            folder / f"smjer_{name}.csv",
            [
                "atom",
                "element",
                "e_x",
                "e_y",
                "e_z",
                "u_po_Q_x",
                "u_po_Q_y",
                "u_po_Q_z",
                "u_x_rel",
                "u_y_rel",
                "u_z_rel",
            ],
            [
                (j, symbol, *e, *u_q, *u_a)
                for j, (symbol, e, u_q, u_a) in enumerate(
                    zip(g4["symbols"], vec, unit_q, item["u_rel"]), 1
                )
            ],
        )
    qe.save_json(folder / "smjerovi_M.json", info)
    return info


def frozen_curve(c, t, g4, name, info, f0):
    """Jedan presjek slobodne energije; A je potpisani najveći pomak atoma u Å."""
    results = tdir(c, t) / "pes" / name
    u = np.asarray(info["u_rel"])
    meff = float(np.sum(mod.masses(g4)[:, None] * u**2))
    # Pretvorba THz² puta amu u zakrivljenost izraženu u eV/Å².
    factor = atomic_mass * (2 * np.pi * 1e12) ** 2 * 1e-20 / elementary_charge
    phonon = factor * meff * info["lambda_THz2"]
    energies = {}
    for j, a in enumerate(c["amplitudes_A"]):
        if a == 0:
            energies[a] = f0
            continue
        d = results / f"pomak_{j:02d}"
        moved = qe.displaced(g4, a * u)
        inp = qe.pw_input(moved, c, t, d, energy_only=True)
        energies[a] = qe.energy(qe.step(d, "scf", "pw", inp, kind="energy_only"))
        release_scratch(d, c)
    mod.csv_out(
        results / "potencijal.csv",
        [
            "A_A",
            "Q_sqrt_amu_A",
            "F_eV_2x2",
            "deltaF_meV_2x2",
            "deltaF_meV_TiSe2",
            "harmonijski_meV_2x2",
        ],
        [
            (a, a * np.sqrt(meff), f, 1000 * (f - f0), 250 * (f - f0), 500 * phonon * a**2)
            for a, f in energies.items()
        ],
    )
    samples = list(energies)
    center = samples.index(0)
    a_minus, a_plus = samples[center - 1], samples[center + 1]
    slope_left = (energies[a_minus] - f0) / a_minus
    slope_right = (energies[a_plus] - f0) / a_plus
    direct = 2 * (slope_right - slope_left) / (a_plus - a_minus)
    minimum = min(samples, key=energies.get)
    minima = [
        {"A_A": a, "deltaF_meV_2x2": 1000 * (energies[a] - f0)}
        for left, a, right in zip(samples, samples[1:], samples[2:])
        if energies[a] < min(energies[left], energies[right])
    ]
    qe.save_json(
        results / "zakrivljenost.json",
        {
            "smjer": name,
            "h_A": a_plus,
            "meff_amu": meff,
            "lambda_THz2": info["lambda_THz2"],
            "rezidual_THz2": info["rezidual_THz2"],
            "PES_eV_A2": direct,
            "DFPT_eV_A2": phonon,
            "najveca_asimetrija_meV": 1000
            * max(abs(energies[a] - energies[b]) for a, b in zip(samples, reversed(samples))),
            "A_najnize_tocke_A": minimum,
            "dubina_na_mrezi_meV_2x2": 1000 * (f0 - energies[minimum]),
            "minimum_na_rubu_intervala": minimum in (samples[0], samples[-1]),
            "unutarnji_minimumi_na_mrezi": minima,
        },
    )


def pes(c, t, g):
    """Dva jednodimenzionalna puta kroz iste položaje pri A=0."""
    info = directions(c, t, g)
    g4 = qe.supercell(g)
    f0 = read_json(tdir(c, t) / "gamma_2x2_ref" / "zavrseno.json")["F_eV"]
    paths = ["M1", "3Q"]
    if info["dodatni_sektor"]:
        paths.append(info["dodatni_sektor"])
    for name in paths:
        frozen_curve(c, t, g4, name, info["smjerovi"][name], f0)


def lock(c, t):
    d = tdir(c, t)
    d.mkdir(parents=True, exist_ok=True)
    return FileLock(d / ".izvodenje.lock")


@contextmanager
def periodic_status(c, root):
    """Ispisuj status i dok glavna nit čeka završetak jednog QE koraka."""
    from izvjestaj import NAMES, status

    log_path = root / "rezultati" / "status.log"
    started = monotonic()
    stop = Event()

    def snapshot(label):
        stamp = datetime.now().astimezone().isoformat(sep=" ", timespec="seconds")
        elapsed = timedelta(seconds=int(monotonic() - started))
        heading = f"[{stamp}] {label} | PID {os.getpid()} | proteklo {elapsed}"
        try:
            disk = shutil.disk_usage(root)
            heading += f" | slobodno {disk.free / 2**30:.2f} od {disk.total / 2**30:.2f} GiB"
        except OSError:
            pass
        try:
            rows = status(c, output=log_path.parent)
            table = "\n".join(
                " ".join(f"{value!s:>12}" for value in row) for row in [["T_K", *NAMES], *rows]
            )
        except OSError as exc:
            table = f"Status trenutačno nije dostupan ({exc})."
        block = f"\n{heading}\n{table}\n"
        try:
            with log_path.open("a", encoding="utf-8") as logfile:
                logfile.write(block)
        except OSError as exc:
            block += f"Upis u {log_path} nije uspio ({exc}).\n"
        try:
            print(block, end="", flush=True)
        except (OSError, ValueError):
            # Odspojen terminal može odbiti ispis. Status je već poslan u datoteku.
            pass

    def repeat():
        while not stop.wait(STATUS_INTERVAL_SECONDS):
            snapshot("Periodični status")

    snapshot(f"Početak praćenja | log {log_path}")
    timer = Thread(target=repeat, name="dft-status", daemon=True)
    timer.start()
    try:
        yield
    finally:
        stop.set()
        timer.join()
        snapshot("Kraj praćenja")


def run(c, t, part, workers):
    with lock(c, t):
        g = relax(c, t, qe.primitive(c), tdir(c, t) / "relax_1x1")
        if part in ("m", "mesh"):
            primitive_part(c, t, g, part)
        elif part == "ref":
            directions(c, t, g)
        elif part == "pes":
            pes(c, t, g)
        elif part == "rel":
            relaxed(c, t, g)
        elif part == "all":
            errors = []

            def attempt(label, work):
                diagnostic = tdir(c, t) / f"greska_{label}.json"
                try:
                    work()
                except (OSError, ValueError, RuntimeError) as exc:
                    qe.save_json(diagnostic, {"faza": label, "poruka": str(exc)})
                    errors.append(f"{label}: {exc}")
                else:
                    diagnostic.unlink(missing_ok=True)

            def normal():
                attempt("mesh", lambda: primitive_part(c, t, g, "mesh"))
                attempt("pes", lambda: pes(c, t, g))

            def distorted_branch():
                attempt("rel", lambda: relaxed(c, t, g))

            if workers == 2:
                with ThreadPoolExecutor(max_workers=2) as pool:
                    fut = pool.submit(distorted_branch)
                    normal()
                    fut.result()
            else:
                normal()
                distorted_branch()
            if errors:
                raise RuntimeError("\n".join(errors))


def prepare(c):
    # Čitljivi predlošci s početnim koordinatama. Run gradi stvarne ulaze nakon ovisnih računa.
    root = qe.run_root(c)
    qe.save_json(root / "postavke.json", c)
    mod.csv_out(
        root / "temperature.csv",
        ["T_el_K", "degauss_Ry", "kBT_eV"],
        [(t, t * qe.KB_EV / qe.RY_EV, t * qe.KB_EV) for t in c["temperatures_K"]],
    )
    g = qe.primitive(c)
    g4 = qe.supercell(g)
    for t in c["temperatures_K"]:
        d = tdir(c, t) / "primjeri_ulaza"
        qe.atomic_text(d / "relax_1x1.in", qe.pw_input(g, c, t, d, relax=True))
        qe.atomic_text(d / "gamma_2x2_scf.in", qe.pw_input(g4, c, t, d))
        qe.atomic_text(
            d / "relax_2x2_random.in", qe.pw_input(qe.distorted(g4, c), c, t, d, relax=True)
        )
        for name, text in [
            ("ph_mesh.in", qe.ph_input(c, "mesh")),
            ("q2r.in", qe.Q2R),
            ("matdyn.in", qe.matdyn_input(c)),
            ("ph_gamma.in", qe.ph_input(c, "gamma")),
            ("dynmat.in", qe.DYNMAT),
            ("ph_M.in", qe.ph_input(c, "pointM")),
        ]:
            qe.atomic_text(d / name, text)
    return root


def run_order(temperatures, part, workers):
    if part == "all" and qe.low_disk():
        # Uz jednog radnika greška bilo koje faze može zaustaviti seriju prije sljedeće faze.
        stages = ("m", "mesh", "pes", "rel") if workers == 1 else ("m", "all")
        return ((stage, t) for t in temperatures for stage in stages)
    stages = ("m", "all") if part == "all" else (part,)
    return ((stage, t) for stage in stages for t in temperatures)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "run", "report", "status", "watch"])
    parser.add_argument("--config", type=Path, default=qe.ROOT / "postavke.json")
    parser.add_argument("--temperatures", nargs="+", type=int, help="Npr. 100 200 400")
    parser.add_argument(
        "--part", choices=["all", "relax1", "m", "mesh", "ref", "pes", "rel"], default="all"
    )
    parser.add_argument(
        "--workers", type=int, choices=[1, 2], default=1, help="Neovisne grane za all"
    )
    args = parser.parse_args()
    c = settings(args.config)
    ts = sorted(set(args.temperatures or c["temperatures_K"]))
    if any(t not in c["temperatures_K"] for t in ts):
        parser.error("Odaberi temperaturu iz postavki.")
    if args.action == "status":
        from izvjestaj import status

        status(c)
        return
    if args.action == "report":
        from izvjestaj import report

        report(c)
        return
    if args.action == "watch":
        from izvjestaj import watch

        watch(c)
        return
    qe.verify_pseudo()
    root = prepare(c)
    print(root, flush=True)
    if args.action == "run":
        from izvjestaj import report

        guard = (
            FileLock(root / ".cijela_serija.lock", timeout=0)
            if args.part == "all"
            else nullcontext()
        )
        problems = []
        try:
            with guard, periodic_status(c, root):
                if args.part == "all":
                    qe.save_json(root / "problemi.json", [])
                for part, t in run_order(ts, args.part, args.workers):
                    try:
                        run(c, t, part, args.workers)
                    except (OSError, ValueError, RuntimeError) as exc:
                        problems.append({"T_K": t, "faza": part, "poruka": str(exc)})
                        print(f"{t} K / {part}: {exc}", flush=True)
                        qe.save_json(root / "problemi.json", problems)
                        if qe.low_disk():
                            break
                    report(c)
        except Timeout:
            raise SystemExit("Cijela serija za ove postavke već je pokrenuta.")
        qe.save_json(root / "problemi.json", problems)
        if problems:
            raise SystemExit(
                "Neke faze zahtijevaju pregled. Vidi problemi.json i odgovarajući QE ispis."
            )


if __name__ == "__main__":
    main()

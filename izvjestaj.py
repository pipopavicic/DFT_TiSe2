"""Tablice i slike nastaju isključivo iz završenih računa."""

import re
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from filelock import FileLock, Timeout

import qe
import modovi as mod
from studija import read_json, tdir


STAGES = [
    ("relax_1x1", "geometrija.json", "relax.out"),
    ("fononi_1x1", "M_direktno.json", "M.out"),
    ("fononi_1x1", "minimum.json", "mesh.out"),
    ("gamma_2x2_ref", "zavrseno.json", "scf.out"),
    ("pes/M1", "zakrivljenost.json", "pomak_*/scf.out"),
    ("pes/3Q", "zakrivljenost.json", "pomak_*/scf.out"),
    ("relax_2x2", "geometrija.json", "relax.out"),
    ("gamma_2x2_rel", "zavrseno.json", "scf.out"),
]
NAMES = [
    "relax_1x1",
    "M_izravno",
    "disperzija",
    "gamma_ref",
    "PES_M1",
    "PES_3Q",
    "relax_2x2",
    "gamma_rel",
]


def status(c, output=None):
    rows = []
    for t in c["temperatures_K"]:
        p = tdir(c, t)
        states = []
        for folder, marker, log in STAGES:
            d = p / folder
            states.append(
                "gotovo"
                if (d / marker).exists()
                else ("započeto" if d.is_dir() and any(d.glob(log)) else "čeka")
            )
        rows.append([t, *states])
    if output:
        mod.csv_out(output / "status.csv", ["T_K", *NAMES], rows)
    else:
        for line in [["T_K", *NAMES], *rows]:
            print(" ".join(f"{str(x):>12}" for x in line))
    return rows


def delta_positions(g, ref):
    cell = np.array(ref["cell_A"])
    frac = (np.array(g["positions_A"]) - np.array(ref["positions_A"])) @ np.linalg.inv(cell)
    delta = (frac - np.rint(frac)) @ cell
    return delta - np.average(delta, axis=0, weights=mod.masses(g))


def align_copies(delta, previous):
    base = np.arange(12).reshape(2, 2, 3)
    candidates = [np.roll(base, (i, j), axis=(0, 1)).ravel() for i, j in qe.COPIES]
    return min(candidates, key=lambda p: np.linalg.norm(delta[p] - previous))


def rel_branches(c, result):
    records = []
    distortion_rows = []
    prev_vec = prev_delta = None
    for t in c["temperatures_K"]:
        p = tdir(c, t)
        rel = p / "gamma_2x2_rel"
        if not (rel / "zavrseno.json").exists():
            continue
        g = read_json(p / "relax_2x2" / "geometrija.json")
        high = qe.supercell(read_json(p / "relax_1x1" / "geometrija.json"))
        delta = delta_positions(g, high)
        with np.load(rel / "gamma.npz") as z:
            f, e = mod.optical(z["nu_THz"], z["e_massweighted"], g)
        if prev_vec is None:
            overlap = np.ones(len(f))
        else:
            perm = align_copies(delta, prev_delta)
            delta = delta[perm]
            e = e.reshape(12, 3, -1)[perm].reshape(36, -1)
            f, e, overlap = mod.match(prev_vec, f, e)
        prev_vec, prev_delta = e, delta
        order = (delta * np.sqrt(mod.masses(g))[:, None]).ravel()
        norm = np.linalg.norm(order)
        projection = (e.T @ order / norm) ** 2 if norm > 0.01 else np.full(len(f), np.nan)
        # Fiksni izbor ishodišta i relativnih predznaka u praćenom referentnom tripletu.
        # Pomak i modovi relaksirane strukture već su usklađeni istom permutacijom kopija.
        triplet = np.full(len(f), np.nan)
        coherent = np.full(len(f), np.nan)
        ref_file = p / "gamma_2x2_ref" / "smjerovi_M.json"
        if ref_file.is_file():
            ref = read_json(ref_file)["smjerovi"]
            arms = np.column_stack([ref[name]["e_massweighted"] for name in mod.M_ARMS])
            triple = np.asarray(ref["3Q"]["e_massweighted"])
            triplet = np.sum((e.T @ arms) ** 2, axis=1)
            coherent = (e.T @ triple) ** 2
            coeff = order @ arms
            total = float(coeff @ coeff / norm**2) if norm > 0.01 else None
            distortion_rows.append((t, norm, *coeff, total))
        gap = np.abs(f[:, None] - f)
        np.fill_diagonal(gap, np.inf)
        near = gap.min(axis=1) < c["instability_threshold_THz"]
        records.append((t, f, overlap, projection, near, triplet, coherent))
    if not records:
        return
    mod.csv_out(
        result / "distorzija_M.csv",
        ["T_K", "norma_sqrt_amu_A", "Q_M1", "Q_M2", "Q_M3", "udio_u_M_tripletu"],
        distortion_rows,
    )
    mod.csv_out(
        result / "grane_relaksirane.csv",
        [
            "T_K",
            "grana",
            "nu_THz",
            "nu2_THz2",
            "preklop_s_prethodnom",
            "udio_u_distorziji",
            "bliska_degeneracija",
            "udio_u_M_tripletu",
            "udio_uz_3Q",
        ],
        [
            (
                t,
                i + 1,
                fi,
                fi * abs(fi),
                ("" if t == records[0][0] else ov[i]),
                (proj[i] if np.isfinite(proj[i]) else ""),
                bool(near[i]),
                triplet[i] if np.isfinite(triplet[i]) else "",
                coherent[i] if np.isfinite(coherent[i]) else "",
            )
            for t, f, ov, proj, near, triplet, coherent in records
            for i, fi in enumerate(f)
        ],
    )
    if len(records) < 2:
        return
    ts = np.array([row[0] for row in records])
    fs = np.array([row[1] for row in records])
    ranges = np.ptp(fs, axis=0)
    win = int(np.argmax(ranges))
    ov = [float(r[2][win]) for r in records[1:]]
    near = any(bool(r[4][win]) for r in records)
    qe.save_json(
        result / "najveci_pomak.json",
        {
            "grana": win + 1,
            "raspon_THz": float(ranges[win]),
            "T_K": ts.tolist(),
            "nu_THz": fs[:, win].tolist(),
            "sve_temperature": len(records) == len(c["temperatures_K"]),
            "najmanji_preklop": min(ov),
            "bliska_degeneracija": near,
            "oznaka_za_provjeru": bool(min(ov) < 0.65 or near),
        },
    )
    fig, ax = plt.subplots(figsize=(7.6, 4.6), layout="constrained")
    ax.plot(ts, fs, color="#a3b2ba", lw=0.7)
    ax.plot(ts, fs[:, win], "o-", color="#a02a46", lw=2, label=f"Grana {win + 1} · najveći raspon")
    ax.axhline(0, color="k", lw=0.8)
    ax.legend()
    ax.set(
        xlabel="Elektronska temperatura (K)",
        ylabel="Predznačena frekvencija (THz)",
        title="Fononi Γ relaksirane superćelije 2×2",
    )
    fig.savefig(result / "relaksirani_modovi.png", dpi=165)
    plt.close(fig)


def mode_illustration(c, dest):
    """Shematski prikaz faze iz geometrije, dostupan prije skupog DFPT računa."""
    g = qe.primitive(c)
    cell = np.array(g["cell_A"])[:2, :2]
    fig, ax = plt.subplots(figsize=(7.1, 4.5), layout="constrained")
    for i, j in qe.COPIES:
        origin = np.array([i, j]) @ cell
        corners = np.array([[0, 0], [1, 0], [1, 1], [0, 1]]) @ cell + origin
        color = "#d6eafa" if i == 0 else "#fbe1cd"
        ax.fill(*corners.T, facecolor=color, edgecolor="#3b4c59", lw=1.1)
        center = np.array([0.5, 0.5]) @ cell + origin
        ax.text(*center, "+" if i == 0 else "−", ha="center", va="center", fontsize=24)
        for sym, r in zip(g["symbols"], g["positions_A"]):
            xy = origin + np.array(r[:2])
            ax.scatter(
                *xy,
                s=48 if sym == "Ti" else 25,
                color="#315275" if sym == "Ti" else "#b75b2e",
                zorder=3,
            )
    for s, color in [("Ti", "#315275"), ("Se", "#b75b2e")]:
        ax.scatter([], [], s=32, color=color, label=s)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(loc="upper right", frameon=False)
    ax.set_title("Za M₁ faza mijenja predznak nakon pomaka za a₁")
    fig.savefig(dest, dpi=160)
    plt.close(fig)


def report(c):
    root = qe.run_root(c)
    result = root / "rezultati"
    result.mkdir(parents=True, exist_ok=True)
    with FileLock(result / ".izvjestaj.lock"):
        _report(c, result)


def watch(c):
    """Osvježi grafove čim pojedina faza završi, pa izađi zajedno sa serijom."""
    root = qe.run_root(c)
    previous = None
    while True:
        snapshot = tuple(
            ((tdir(c, t) / folder / marker).exists(), any((tdir(c, t) / folder).glob(log)))
            for t in c["temperatures_K"]
            for folder, marker, log in STAGES
        )
        if snapshot != previous:
            report(c)
            previous = snapshot
        try:
            with FileLock(root / ".cijela_serija.lock", timeout=0):
                report(c)
                return
        except Timeout:
            time.sleep(45)


def in_plane_stress(directory):
    """Tenzor naprezanja ploče u konvenciji ASE, s pozitivnim vlakom, u eV/Å²."""
    lines = (directory / "relax.out").read_text(errors="replace").splitlines()
    headers = [i for i, line in enumerate(lines) if "stress" in line and "(Ry/bohr**3)" in line]
    if not headers:
        raise ValueError(f"Nedostaje završni tenzor naprezanja u {directory}.")
    start = headers[-1] + 1
    pressure = np.array([mod.numbers(line)[:3] for line in lines[start : start + 3]])
    cell = np.array(read_json(directory / "geometrija.json")["cell_A"])
    height = abs(np.linalg.det(cell)) / np.linalg.norm(np.cross(cell[0], cell[1]))
    return -pressure[:2, :2] * (qe.RY_EV / qe.BOHR_A**3) * height


def pes_report(curves, dest):
    """Usporedi putove po istoj fizičkoj amplitudi i po istoj superćeliji."""
    rows = []
    axis_labels = {
        "xlabel": "Potpisani najveći atomski pomak A (Å)",
        "ylabel": "F(A) − F(0) (meV po ćeliji 2×2)",
    }
    for t, paths in curves.items():
        for name, (data, stats) in paths.items():
            minima = stats["unutarnji_minimumi_na_mrezi"]
            rows.append(
                (
                    t,
                    name,
                    stats["A_najnize_tocke_A"],
                    stats["dubina_na_mrezi_meV_2x2"],
                    stats["minimum_na_rubu_intervala"],
                    stats["PES_eV_A2"],
                    stats["DFPT_eV_A2"],
                    stats["najveca_asimetrija_meV"],
                    " ".join(f"{point['A_A']:g}" for point in minima),
                )
            )
        if not {"M1", "3Q"}.issubset(paths):
            continue
        fig, ax = plt.subplots(figsize=(7.5, 4.8), layout="constrained")
        for name, label, marker, color in [
            ("M1", "Samo M₁", "o", "#225989"),
            ("3Q", "M₁ + M₂ + M₃", "s", "#993e52"),
        ]:
            data = paths[name][0]
            ax.plot(data[:, 0], data[:, 3], marker=marker, ms=4, color=color, label=label)
        ax.axhline(0, color="0.45", lw=0.7)
        ax.set(**axis_labels, title=f"TiSe₂ pri {t} K · dva presjeka slobodne energije")
        ax.legend()
        fig.savefig(dest / f"PES_M1_3Q_{t:04d}K.png", dpi=165)
        plt.close(fig)
    mod.csv_out(
        dest / "PES_usporedba.csv",
        [
            "T_K",
            "smjer",
            "A_najnize_tocke_A",
            "dubina_na_mrezi_meV_2x2",
            "minimum_na_rubu",
            "PES_eV_A2",
            "DFPT_eV_A2",
            "asimetrija_meV_2x2",
            "A_unutarnjih_minimuma_A",
        ],
        rows,
    )
    for name, filename in [("M1", "PES_prema_T.png"), ("3Q", "PES_3Q_prema_T.png")]:
        series = [(t, paths[name][0]) for t, paths in curves.items() if name in paths]
        if series:
            fig, ax = plt.subplots(figsize=(7.5, 4.8), layout="constrained")
            for t, data in series:
                ax.plot(data[:, 0], data[:, 3], "o-", ms=3, label=f"{t} K")
            ax.axhline(0, color="0.45", lw=0.7)
            ax.set(**axis_labels, title=f"Presjek {name} pri različitim elektronskim temperaturama")
            ax.legend(ncol=2)
            fig.savefig(dest / filename, dpi=165)
            plt.close(fig)


def mode_arrows(g, directions, t, dest):
    """Tlocrt DFT pomaka u referentnoj superćeliji, s istim mjerilom za oba smjera."""
    cell = np.asarray(g["cell_A"])
    frac = np.asarray(g["positions_A"]) @ np.linalg.inv(cell)
    frac[:, :2] = (frac[:, :2] + 1e-7) % 1 - 1e-7
    pos = frac @ cell
    xy = pos[:, :2]
    ti = np.asarray(g["symbols"]) == "Ti"
    upper = pos[:, 2] > pos[ti, 2].mean()
    outline = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]) @ cell[:2, :2]
    for name in ("M1", "3Q"):
        arrows = 0.8 * np.asarray(directions[name]["u_rel"])[:, :2]
        fig, ax = plt.subplots(figsize=(6.8, 5.2), layout="constrained")
        ax.plot(*outline.T, lw=1, color="0.6")
        for mask, marker, color, label in [
            (ti, "o", "#3674a5", "Ti"),
            (~ti & upper, "^", "#daa353", "Se iznad Ti"),
            (~ti & ~upper, "v", "#b66b30", "Se ispod Ti"),
        ]:
            ax.scatter(*xy[mask].T, marker=marker, color=color, s=55, label=label)
        ax.quiver(*xy.T, *arrows.T, angles="xy", scale_units="xy", scale=1, color="#28303a")
        ax.update_datalim(xy + arrows)
        ax.margins(0.13)
        ax.set_aspect("equal")
        ax.set(
            xlabel="x (Å)",
            ylabel="y (Å)",
            title=f"TiSe₂ · {name} pri {t} K\nProjekcija pomaka u ravninu, uvećane strelice",
        )
        ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        fig.savefig(dest / f"pomaci_{name}_{t:04d}K.png", dpi=165)
        plt.close(fig)


def _report(c, result):
    state = status(c, result)
    summary = []
    mpoints = []
    triplet_rows = []
    curves = {}
    geometry = []
    stress = []
    for t in c["temperatures_K"]:
        p = tdir(c, t)
        for cell_name in ["relax_1x1", "relax_2x2"]:
            directory = p / cell_name
            if (directory / "geometrija.json").exists() and (directory / "relax.ok.json").exists():
                sigma = in_plane_stress(directory)
                stress.append((t, cell_name, sigma[0, 0], sigma[0, 1], sigma[1, 1]))
        r = p / "relax_1x1"
        if (r / "geometrija.json").exists() and (r / "relax.ok.json").exists():
            g = read_json(r / "geometrija.json")
            z = np.array(g["positions_A"])[:, 2]
            forces = (
                (r / "relax.out")
                .read_text(errors="replace")
                .rsplit("Forces acting on atoms", 1)[-1]
            )
            vectors = re.findall(r"atom\s+\d+\s+type\s+\d+\s+force\s*=\s*([^\n]+)", forces)
            fmax = (
                max(np.linalg.norm(mod.numbers(row)) for row in vectors) * qe.RY_EV / qe.BOHR_A
                if vectors
                else None
            )
            geometry.append((t, qe.energy(r / "relax.out"), abs(z[1] - z[2]) / 2, fmax))
        d = p / "fononi_1x1"
        minimum = read_json(d / "minimum.json") if (d / "minimum.json").exists() else {}
        m = minimum.get("nu_M_THz")
        direct = (
            read_json(d / "M_direktno.json")["nu_min_THz"]
            if (d / "M_direktno.json").exists()
            else None
        )
        if minimum:
            data = np.loadtxt(d / "disperzija.csv", delimiter=",", skiprows=1, ndmin=2)
            xx, yy = data[:, 0], data[:, 4:]
            fig, ax = plt.subplots(figsize=(7.5, 4.7), layout="constrained")
            ax.plot(xx, yy, color="#244f74", lw=1.2)
            ax.plot(
                xx,
                np.ma.masked_where(yy >= -c["instability_threshold_THz"], yy),
                color="#a32948",
                lw=1.7,
            )
            ax.fill_between(
                xx, yy.min(axis=1), 0, where=yy.min(axis=1) < 0, color="#eab9bd", interpolate=True
            )
            ax.axhline(0, color="k", lw=0.8)
            ticks = xx[:: c["path_steps_per_segment"]]
            ax.set_xticks(ticks, ["Γ", "M", "K", "Γ"])
            ax.grid(axis="x", alpha=0.3)
            ax.set(
                ylabel="Predznačena frekvencija (THz)",
                title=f"TiSe₂ · DFT-PBE · {t} K\nk = {c['k_mesh']}×{c['k_mesh']}×1, q = {c['q_mesh']}×{c['q_mesh']}×1",
            )
            fig.savefig(result / f"disperzija_{t:04d}K.png", dpi=165)
            plt.close(fig)
        ref = p / "gamma_2x2_ref"
        folded = None
        if (ref / "zavrseno.json").exists():
            with np.load(ref / "gamma.npz") as z:
                folded = float(mod.sector_modes(z["D_signed_THz2"])["M1"][0][0])
        if (ref / "smjerovi_M.json").is_file():
            directions = read_json(ref / "smjerovi_M.json")["smjerovi"]
            triple = directions["3Q"]
            triplet_rows.append(
                (
                    t,
                    *(directions[k]["nu_THz"] for k in mod.M_ARMS),
                    triple["lambda_THz2"],
                    triple["raspon_lambda_M_THz2"],
                    triple["rezidual_THz2"],
                )
            )
            high = qe.supercell(read_json(p / "relax_1x1" / "geometrija.json"))
            mode_arrows(high, directions, t, result)
        if m is not None or direct is not None or folded is not None:
            mpoints.append((t, m, direct, folded))
        rp = p / "gamma_2x2_rel"
        diff = None
        relmin = None
        if (rp / "zavrseno.json").exists():
            with np.load(rp / "gamma.npz") as z:
                relmin = float(z["nu_THz"].min())
            if (ref / "zavrseno.json").exists():
                diff = 1e3 * (
                    read_json(rp / "zavrseno.json")["F_eV"]
                    - read_json(ref / "zavrseno.json")["F_eV"]
                )
        for marker in sorted((p / "pes").glob("*/zakrivljenost.json")):
            stats = read_json(marker)
            data = np.loadtxt(
                marker.with_name("potencijal.csv"), delimiter=",", skiprows=1, ndmin=2
            )
            curves.setdefault(t, {})[stats["smjer"]] = (data, stats)
        summary.append(
            (
                t,
                t * qe.KB_EV / qe.RY_EV,
                m,
                direct,
                folded,
                None if folded is None or m is None else folded - m,
                relmin,
                diff,
            )
        )
    mod.csv_out(
        result / "M_triplet_prema_T.csv",
        [
            "T_K",
            "nu_M1_THz",
            "nu_M2_THz",
            "nu_M3_THz",
            "lambda_3Q_THz2",
            "raspon_lambda_M_THz2",
            "rezidual_3Q_THz2",
        ],
        triplet_rows,
    )
    mod.csv_out(
        result / "naprezanje_2d.csv",
        ["T_K", "struktura", "s_xx_eV_A2", "s_xy_eV_A2", "s_yy_eV_A2"],
        stress,
    )
    mod.csv_out(
        result / "relaksacija_1x1.csv",
        ["T_K", "F_eV_1x1", "visina_Se_iznad_Ti_A", "Fmax_eV_A"],
        geometry,
    )
    mod.csv_out(
        result / "usporedba.csv",
        [
            "T_K",
            "degauss_Ry",
            "M_putanja_THz",
            "M_direktno_THz",
            "M1_u_Gamma_THz",
            "razlika_2x2_minus_1x1_THz",
            "min_relaksirana_Gamma_THz",
            "deltaF_rel_minus_ref_meV_2x2",
        ],
        summary,
    )
    if mpoints:
        fig, ax = plt.subplots(figsize=(7.5, 4.5), layout="constrained")
        for j, label, marker in [
            (1, "M iz disperzije", "o-"),
            (2, "Izravno M", "x--"),
            (3, "M₁ u Γ superćelije", "s:"),
        ]:
            vals = [(row[0], row[j]) for row in mpoints if row[j] is not None]
            if vals:
                ax.plot(*zip(*vals), marker, ms=5, label=label)
        ax.axhline(0, color="k", lw=0.8)
        ax.legend()
        ax.set(
            xlabel="Elektronska temperatura (K)",
            ylabel="Predznačena frekvencija (THz)",
            title="Najniža grana u M",
        )
        fig.savefig(result / "M_prema_T.png", dpi=165)
        plt.close(fig)
    pes_report(curves, result)
    rel_branches(c, result)
    lines = [
        "# Izračunati rezultati",
        "",
        f"Elektronska mreža {c['k_mesh']}×{c['k_mesh']}×1, fononska mreža {c['q_mesh']}×{c['q_mesh']}×1.",
        "",
        "| T (K) | Relaksacija 1×1 | M izravno | Disperzija | Γ simetrične 2×2 | PES M1 | PES 3Q | Relaksacija 2×2 | Γ relaksirane 2×2 |",
        "|---:|---|---|---|---|---|---|---|---|",
    ]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in state]
    if geometry:
        lines += [
            "",
            "## Visokosimetrična primitivna ćelija",
            "",
            "| T (K) | Visina Se iznad Ti (Å) | Slobodna energija (eV/ćeliji) | Najveća sila (eV/Å) |",
            "|---:|---:|---:|---:|",
        ]
        for t, f, h, fmax in geometry:
            val = "—" if fmax is None else f"{fmax:.6g}"
            lines.append(f"| {t} | {h:.6f} | {f:.8f} | {val} |")
    if mpoints:
        lines += [
            "",
            "## Najniža predznačena frekvencija u M (THz)",
            "",
            "| T (K) | Iz disperzije | Izravni DFPT | M₁ presavijen u Γ |",
            "|---:|---:|---:|---:|",
        ]
        for t, *freq in mpoints:
            lines.append(
                "| "
                + str(t)
                + " | "
                + " | ".join("—" if f is None else f"{f:.5f}" for f in freq)
                + " |"
            )
    if curves:
        lines += [
            "",
            "## Dubina presjeka na uzorkovanim amplitudama (meV po 2×2)",
            "",
            "| T (K) | M1 | 3Q |",
            "|---:|---:|---:|",
        ]
        for t, paths in curves.items():
            values = []
            for name in ("M1", "3Q"):
                if name in paths:
                    stats = paths[name][1]
                    edge = " (rub)" if stats["minimum_na_rubu_intervala"] else ""
                    values.append(f"{stats['dubina_na_mrezi_meV_2x2']:.3f}{edge}")
                else:
                    values.append("—")
            lines.append(f"| {t} | {values[0]} | {values[1]} |")
    qe.atomic_text(result / "PREGLED.md", "\n".join(lines) + "\n")

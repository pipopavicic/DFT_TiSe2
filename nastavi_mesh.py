"""Nastavak prekinute disperzije iz njezinih spremljenih matrica (Srce, Linux)."""

import argparse
from contextlib import ExitStack
import filecmp
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile


def q_grid(path):
    import numpy as np

    grid, count, *rows = path.read_text().splitlines()
    points = np.array([list(map(float, row.split())) for row in rows if row.strip()])
    if points.shape != (int(count), 3) or not np.isfinite(points).all():
        raise ValueError(f"Nepotpuna q-mreža u {path}")
    return tuple(map(int, grid.split())), points


def complete_matrix(path, q):
    """Provjeri q pri dijagonalizaciji, svih devet frekvencija i kraj zapisa."""
    import numpy as np

    if not path.is_file():
        return False
    text = path.read_text(errors="replace")
    _, header, modes = text.partition("Diagonalizing the dynamical matrix")
    position = re.search(r"q\s*=\s*\(([^)]*)\)", modes)
    if not header or not position or not modes.rstrip().endswith("********"):
        return False
    freq = re.findall(r"^\s*freq\s*\(\s*(\d+)\s*\)\s*=\s*(\S+)\s*\[THz\]", modes, re.M)
    try:
        coords = np.array([float(x) for x in position[1].split()])
        finite = all(np.isfinite(float(nu.lower().replace("d", "e"))) for _, nu in freq)
    except ValueError:
        return False
    return (coords.shape == (3,) and np.allclose(coords, q, rtol=0, atol=1e-6)
            and [int(i) for i, _ in freq] == list(range(1, 10))
            and finite)


def offload_old_ph(directory, scratch, t):
    """Premjesti velike radne datoteke nakon usporedbe izvornika i kopije."""
    old = directory / "tmp" / "_ph0"
    if old.resolve() != old or scratch.is_relative_to(directory):
        raise ValueError("Očekujem izvornu mapu tmp/_ph0 na disku s rezultatima i zaseban scratch.")
    for proc in Path("/proc").glob("[0-9]*"):
        try:
            if (proc / "comm").read_text().strip() in {"pw.x", "ph.x"} and (
                Path(os.readlink(proc / "cwd")).is_relative_to(directory)
            ):
                raise RuntimeError(f"QE još radi u {directory}.")
        except OSError:
            continue
    files = []
    for parent, subdirs, names in os.walk(old):
        subdirs[:] = [s for s in subdirs if not s.endswith(".phsave")]
        for name in names:
            path = Path(parent) / name
            if not path.is_symlink() and path.is_file() and path.stat().st_size >= 1024**2:
                files.append(path)
    if files:
        archive = Path(tempfile.mkdtemp(prefix=f"tise2-{t}K-prekid-", dir=scratch))
        print(f"Kopija starih radnih datoteka: {archive}", flush=True)
        for src in sorted(files, key=lambda p: p.stat().st_size, reverse=True):
            dest = archive / src.relative_to(old)
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                with src.open("rb") as original, dest.open("xb") as copy:
                    shutil.copyfileobj(original, copy, length=1024**2)
                    copy.flush()
                    os.fsync(copy.fileno())
                shutil.copystat(src, dest)
                if not filecmp.cmp(src, dest, shallow=False):
                    raise OSError(f"Kopija se razlikuje. Izvornik je sačuvan na {src}.")
            except (OSError, KeyboardInterrupt):
                dest.unlink(missing_ok=True)
                raise
            src.unlink()
    print(f"Slobodno uz rezultate: {shutil.disk_usage(directory).free / 2**30:.2f} GiB; "
          f"na radnom disku: {shutil.disk_usage(scratch).free / 2**30:.2f} GiB.", flush=True)


def resumed_input(original, iq, work):
    # Fildyn ostaje mesh.dyn. Svaki ph.x radi u vlastitoj mapi pa ne dira stare matrice.
    quoted = "'" + work.as_posix().replace("'", "''") + "'"
    text, paths = re.subn(r"""(?i)(\boutdir\s*=\s*)['"]\./tmp/?['"]""",
                          lambda match: match[1] + quoted, original)
    if paths != 1 or not re.search(r"""(?i)\bfildyn\s*=\s*['"]mesh\.dyn['"]""", text):
        raise ValueError("Očekujem izvorni mesh.in projekta s outdir='./tmp' i fildyn='mesh.dyn'.")
    if re.search(r"(?i)\b(start_q|last_q|recover)\s*=", text):
        raise ValueError("Izvorni mesh.in već ima odabir točaka ili postavku recover.")
    text, count = re.subn(r"(?i)&INPUTPH\b", f"&INPUTPH\n start_q={iq}, last_q={iq}, recover=.false.,", text)
    if count != 1:
        raise ValueError("Nedostaje jedinstven blok &INPUTPH u mesh.in.")
    return text


def install_file(source, dest):
    """Prethodni izlaz ostaje čitljiv sve do dovršenog kopiranja."""
    with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as f:
        temp = Path(f.name)
    try:
        shutil.copyfile(source, temp)
        temp.replace(dest)
    finally:
        temp.unlink(missing_ok=True)


def point_confirmed(stage, iq, grid, points):
    import numpy as np

    if not (stage / "mesh.ok.json").is_file():
        return False
    try:
        new_grid, new_points = q_grid(stage / "mesh.dyn0")
        return (new_grid == grid and new_points.shape == points.shape
                and np.allclose(new_points, points, rtol=0, atol=1e-8)
                and complete_matrix(stage / f"mesh.dyn{iq}", points[iq - 1]))
    except (FileNotFoundError, ValueError):
        return False


def recover(directory, c, t, scratch, points, original, qe, mod):
    saved = directory / "tmp/tise2.save"
    geometry = json.loads((directory.parent / "relax_1x1/geometrija.json").read_text(encoding="utf-8-sig"))
    steps = directory / "nastavak_q"
    steps.mkdir(exist_ok=True)
    grid, _ = q_grid(directory / "mesh.dyn0")
    for iq, q in enumerate(points, 1):
        target = directory / f"mesh.dyn{iq}"
        if complete_matrix(target, q):
            continue
        stage = steps / f"q{iq:02}"
        stage.mkdir(exist_ok=True)
        if target.exists() and target.stat().st_size:
            # Čuvamo i nepotpuni izvorni zapis prije zamjene.
            if not (stage / "izvorna_nepotpuna.dyn").exists():
                install_file(target, stage / "izvorna_nepotpuna.dyn")
        if not point_confirmed(stage, iq, grid, points):
            print(f"Računam q {iq}/{len(points)}. Izlaz: {stage / 'mesh.out'}", flush=True)
            with tempfile.TemporaryDirectory(prefix=f"tise2-{t}K-q{iq}-", dir=scratch) as temp:
                work = Path(temp)
                shutil.copytree(saved, work / saved.name)
                for name, _ in qe.UPF.values():
                    pseudo = work / saved.name / name
                    if not pseudo.is_file():
                        shutil.copy2(qe.ROOT / "pseudo" / name, pseudo)
                for wfc in saved.parent.glob("tise2.wfc*"):
                    if wfc.is_file():
                        shutil.copy2(wfc, work / wfc.name)
                qe.step(stage, "mesh", "ph", resumed_input(original, iq, work),
                        ["mesh.dyn0", target.name])
            if not point_confirmed(stage, iq, grid, points):
                raise ValueError(f"Mreža ili matrica za q {iq} nije potvrđena. Sačuvan je prethodni izlaz.")
        install_file(stage / target.name, target)
        print(f"Trajno spremljeno: {target.name}", flush=True)

    # Pokreni interpolaciju nad točno ovim skupom matrica, u čistoj mapi.
    with tempfile.TemporaryDirectory(prefix="interpolacija-", dir=steps) as post:
        post = Path(post)
        for i in range(len(points) + 1):
            name = f"mesh.dyn{i}"
            shutil.copy2(directory / name, post / name)
        qe.step(post, "q2r", "q2r", qe.Q2R, ["mesh.fc"])
        qe.step(post, "matdyn", "matdyn", qe.matdyn_input(c), ["path.freq"])
        minimum = mod.dispersion(post, c, geometry)
        for src in post.iterdir():
            if src.is_file() and not src.name.startswith("mesh.dyn"):
                install_file(src, directory / src.name)
        qe.save_json(directory / "minimum.json", minimum)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("racun", type=Path, help="Mapa računa s datotekom postavke.json")
    parser.add_argument("temperatura", type=int)
    parser.add_argument("--scratch", type=Path, default=Path("/tmp"))
    parser.add_argument("--od-q", type=int, default=1,
                        help="Prva točka koju se smije računati. Prethodne moraju biti potpune.")
    parser.add_argument("--premjesti", action="store_true",
                        help="Samo oslobodi trajni disk premještanjem provjerenih kopija radnih datoteka.")
    args = parser.parse_args()
    root, scratch, t = args.racun.resolve(), args.scratch.resolve(strict=True), args.temperatura
    project = root.parent.parent
    sys.path.insert(0, str(project))
    import modovi
    import qe
    import studija
    from filelock import FileLock

    c = studija.settings(root / "postavke.json")
    if qe.run_root(c).resolve() != root or t not in c["temperatures_K"]:
        parser.error("Mapa i temperatura moraju pripadati odabranoj studiji.")
    directory = root / f"T{t:04d}K" / "fononi_1x1"
    if not scratch.is_dir() or scratch.stat().st_dev == directory.stat().st_dev:
        parser.error("Odaberi radnu mapu na drugom datotečnom sustavu.")
    with ExitStack() as stack:
        series_lock = root / ".cijela_serija.lock"
        if series_lock.exists():
            stack.enter_context(FileLock(series_lock, timeout=0))
        stack.enter_context(studija.lock(c, t).acquire(timeout=0))
        grid, points = q_grid(directory / "mesh.dyn0")
        if grid != (c["q_mesh"], c["q_mesh"], 1):
            raise ValueError("Mreža u mesh.dyn0 ne odgovara postavkama studije.")
        original = (directory / "mesh.in").read_text()
        resumed_input(original, 1, scratch)
        good = [i for i, q in enumerate(points, 1) if complete_matrix(directory / f"mesh.dyn{i}", q)]
        print("Potvrđene q-točke:", good, "od", len(points), flush=True)
        if not 1 <= args.od_q <= len(points) or any(i not in good for i in range(1, args.od_q)):
            parser.error(f"Provjeri prvi dopušteni indeks i izvorne matrice prije q {args.od_q}.")
        if len(good) < len(points):
            state = directory / "tmp/tise2.save"
            if not (state / "data-file-schema.xml").is_file() or not (
                any(state.glob("wfc*")) or any(state.parent.glob("tise2.wfc*"))
            ):
                raise ValueError("Nedostaje spremljeni SCF ili njegove valne funkcije u izvornom tmp.")
        if args.premjesti:
            offload_old_ph(directory, scratch, t)
            reader = Path(__file__).resolve().with_name("github_status.py")
            print("Ćelija za nastavak u Python bilježnici:", flush=True)
            print(f'%run "{Path(__file__).resolve()}" "{root}" {t} --od-q {args.od_q} '
                  f'--scratch "{scratch}"', flush=True)
            if reader.is_file():
                print("Ćelija za GitHub u drugoj Python bilježnici nakon pokretanja nastavka:", flush=True)
                print(f'%run "{reader}" "{root}"', flush=True)
            return
        with studija.periodic_status(c, root):
            if not (directory / "minimum.json").exists() or len(good) != len(points):
                recover(directory, c, t, scratch, points, original, qe, modovi)
            from izvjestaj import report
            report(c)
        print(f"Disperzija na {t} K: {directory / 'disperzija.csv'}", flush=True)


if __name__ == "__main__":
    main()

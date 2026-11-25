"""Ulazi za QE i izvođenje jednog koraka s provjerom konvergencije."""

from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
from scipy.constants import physical_constants

ROOT = Path(__file__).resolve().parent
RY_EV = physical_constants["Rydberg constant times hc in eV"][0]
BOHR_A = physical_constants["Bohr radius"][0] * 1e10
KB_EV = physical_constants["Boltzmann constant in eV/K"][0]
MASSES = {"Ti": 47.867, "Se": 78.971}
UPF = {
    "Ti": (
        "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "6f97c81d27b063715dac8214c3d983608c1cc1a2d7f9a4b00eea25c67d1ffd7e",
    ),
    "Se": (
        "Se.pbe-n-kjpaw_psl.1.0.0.UPF",
        "92c94ecdb4217de69da39e4c36e96c178e38da9c7b574ef83df317b16d9e403c",
    ),
}
REAL = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"


def atomic_text(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="\n", dir=path.parent, delete=False
    ) as f:
        f.write(content)
        tmp = Path(f.name)
    tmp.replace(path)


def save_json(path, data):
    atomic_text(path, json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def verify_pseudo():
    for element, (name, expected) in UPF.items():
        if hashlib.sha256((ROOT / "pseudo" / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Pogrešan SHA-256 pseudopotencijala za {element}")


def run_root(config):
    key = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()[:10]
    return ROOT / "racuni" / f"k{config['k_mesh']}_q{config['q_mesh']}_{key}"


def geom(cell, symbols, positions):
    return {
        "cell_A": np.asarray(cell).tolist(),
        "symbols": list(symbols),
        "positions_A": np.asarray(positions).tolist(),
    }


def primitive(c):
    a, h, z = c["a_A"], c["c_A"], c["se_height_A"]
    cell = np.array([[a, 0, 0], [-a / 2, np.sqrt(3) * a / 2, 0], [0, 0, h]])
    frac = [[0, 0, 0.5], [2 / 3, 1 / 3, 0.5 + z / h], [1 / 3, 2 / 3, 0.5 - z / h]]
    return geom(cell, ["Ti", "Se", "Se"], frac @ cell)


COPIES = [(0, 0), (0, 1), (1, 0), (1, 1)]


def supercell(g):
    cell = np.array(g["cell_A"])
    pos = np.array(g["positions_A"])
    copies = [pos + i * cell[0] + j * cell[1] for i, j in COPIES]
    return geom(np.diag([2, 2, 1]) @ cell, g["symbols"] * 4, np.vstack(copies))


def displaced(g, shift):
    return geom(g["cell_A"], g["symbols"], np.array(g["positions_A"]) + shift)


def distorted(g, c):
    masses = np.array([MASSES[s] for s in g["symbols"]])
    shift = np.random.default_rng(c["random_seed"]).normal(size=(len(masses), 3))
    shift -= np.average(shift, axis=0, weights=masses)
    shift *= c["random_rms_A"] / np.sqrt(np.mean(np.sum(shift**2, axis=1)))
    return displaced(g, shift)


def xyz(path, g):
    lattice = " ".join(f"{x:.10f}" for x in np.ravel(g["cell_A"]))
    rows = [
        str(len(g["symbols"])),
        f'Lattice="{lattice}" Properties=species:S:1:pos:R:3 pbc="T T F"',
    ]
    rows += [
        f"{s} " + " ".join(f"{x:.10f}" for x in v) for s, v in zip(g["symbols"], g["positions_A"])
    ]
    atomic_text(path, "\n".join(rows) + "\n")


def pw_input(g, c, t, directory, *, relax=False, energy_only=False):
    n = len(g["symbols"]) // 3
    k = c["k_mesh"] // (2 if n == 4 else 1)
    cell = np.array(g["cell_A"])
    frac = np.linalg.solve(cell.T, np.array(g["positions_A"]).T).T
    pp = Path(os.path.relpath(ROOT / "pseudo", directory)).as_posix().replace("'", "''")
    # Nasumičan početak već razbija simetriju. U nemagnetskom računu zadržavamo
    # vremenski obrat i simetrije koje QE u konkretnoj geometriji stvarno pronađe.
    control = "relax" if relax else "scf"
    io = ("medium" if n == 4 else "low") if not energy_only else ("medium" if n == 4 else "nowf")
    lines = [
        "&CONTROL",
        f" calculation='{control}', prefix='tise2', outdir='./tmp', pseudo_dir='{pp}',",
        f" disk_io='{io}', tprnfor=.true., tstress=.true., nstep=120,",
        f" forc_conv_thr={c['force_threshold_Ry_Bohr']}, etot_conv_thr=1.0d-6,",
        "/",
        "&SYSTEM",
        f" ibrav=0, nat={len(frac)}, ntyp=2, input_dft='PBE',",
        f" ecutwfc={c['ecutwfc_Ry']}, ecutrho={c['ecutrho_Ry']}, nbnd={n * c['bands_per_formula_unit']},",
        f" occupations='smearing', smearing='fd', degauss={KB_EV * t / RY_EV:.13f},",
        "/",
        "&ELECTRONS",
        f" conv_thr={c['conv_thr_Ry']}, electron_maxstep=300, mixing_beta=0.30,",
        "/",
    ]
    if relax:
        lines += ["&IONS", " ion_dynamics='bfgs',", "/"]
    lines += ["ATOMIC_SPECIES"] + [f"{s} {MASSES[s]} {UPF[s][0]}" for s in MASSES]
    lines += ["CELL_PARAMETERS angstrom"] + [" ".join(f"{x:.13f}" for x in row) for row in cell]
    lines += ["ATOMIC_POSITIONS crystal"]
    for s, v in zip(g["symbols"], frac):
        flags = (" 0 0 0" if s == "Ti" else " 0 0 1") if relax and n == 1 else ""
        lines.append(s + " " + " ".join(f"{x:.13f}" for x in v) + flags)
    return "\n".join(lines + ["K_POINTS automatic", f"{k} {k} 1 0 0 0", ""])


def ph_input(c, mode):
    if mode not in ("mesh", "gamma", "pointM"):
        raise ValueError(f"Nepoznata vrsta fononskog računa: {mode}")
    common = [
        "TiSe2 DFPT",
        "&INPUTPH",
        " prefix='tise2', outdir='./tmp',",
        f" tr2_ph={c['phonon_threshold']}, alpha_mix(1)=0.3,",
        " epsil=.false., trans=.true.,",
    ]
    if mode == "mesh":
        q = c["q_mesh"]
        return "\n".join(
            common + [" fildyn='mesh.dyn', ldisp=.true.,", f" nq1={q}, nq2={q}, nq3=1,", "/", ""]
        )
    point = [0.0, 0.0, 0.0] if mode == "gamma" else [0.5, 0.5 / np.sqrt(3), 0.0]
    return "\n".join(
        common
        + [f" fildyn='{mode}.dyn', ldisp=.false.,", "/", " ".join(f"{x:.13f}" for x in point), ""]
    )


def path_points(c):
    vertices = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1 / 3, 1 / 3, 0.0], [0.0, 0.0, 0.0]])
    return np.vstack(
        [
            *(
                np.linspace(a, b, c["path_steps_per_segment"], endpoint=False)
                for a, b in zip(vertices[:-1], vertices[1:])
            ),
            vertices[-1:],
        ]
    )


def matdyn_input(c):
    points = path_points(c)
    lines = [
        "&INPUT",
        " flfrc='mesh.fc', flfrq='path.freq', asr='crystal',",
        " q_in_cryst_coord=.true., q_in_band_form=.false.,",
        "/",
        str(len(points)),
    ]
    return "\n".join(lines + [" ".join(f"{x:.13f}" for x in v) for v in points] + [""])


Q2R = "&INPUT\n fildyn='mesh.dyn', flfrc='mesh.fc', zasr='crystal',\n/\n"
DYNMAT = (
    "&INPUT\n fildyn='gamma.dyn', asr='crystal', fileig='gamma.eig', filout='gamma.modes',\n/\n"
)


def executable(name):
    for directory in [os.environ.get("QE_BIN"), ROOT / "engine"]:
        if directory:
            for suff in (".x", ".exe", ".x.exe"):
                p = Path(directory).resolve() / (name + suff)
                if p.is_file():
                    return str(p)
    for suff in (".x", ".exe", ".x.exe"):
        found = shutil.which(name + suff)
        if found:
            return found
    raise FileNotFoundError(f"{name}.x nije pronađen. Postavi QE_BIN na mapu s programima.")


def low_disk():
    return os.environ.get("TISE2_LOW_DISK") == "1"


def disk_saving_input(contents, program, kind):
    """Promijeni samo način pohrane u ulazima koje generira ova studija."""
    if program == "ph" and "reduce_io" not in contents.lower():
        return contents.replace("&INPUTPH\n", "&INPUTPH\n reduce_io=.true.,\n", 1)
    if program == "pw" and kind in ("phonon_scf", "relax", "energy_only"):
        # Samo SCF iz kojeg slijedi ph.x mora na kraju zapisati valne funkcije.
        io = "low" if kind == "phonon_scf" else "nowf"
        return re.sub(r"(?m)^ disk_io='(?:low|medium|nowf)',", f" disk_io='{io}',", contents, count=1)
    return contents


def step(directory, label, program, contents, needed=(), kind=None):
    """Svaka potvrda vrijedi samo uz odgovarajući ulaz i potrebne izlazne datoteke."""
    directory.mkdir(parents=True, exist_ok=True)
    inp, log, receipt = [directory / (label + s) for s in (".in", ".out", ".ok.json")]
    binary = executable(program)
    binary_signature = f"\n{Path(binary).name}\n{Path(binary).stat().st_size}"
    fingerprint = hashlib.sha256((contents + binary_signature).encode()).hexdigest()
    if receipt.exists() and log.is_file() and all((directory / p).is_file() for p in needed):
        previous = json.loads(receipt.read_text())
        if previous.get("workflow_sha256", previous["input_sha256"]) == fingerprint:
            return log
    receipt.unlink(missing_ok=True)
    actual_input = disk_saving_input(contents, program, kind) if low_disk() else contents
    atomic_text(inp, actual_input)
    cmd = [binary, "-in", inp.name]
    ranks = int(os.environ.get("TISE2_MPI_RANKS", "1"))
    if ranks > 1 and program in ("pw", "ph"):
        cmd = [os.environ.get("MPIEXEC", "mpiexec"), "-n", str(ranks), *cmd]
    env = dict(os.environ)
    threads = env.get("TISE2_THREADS", "2")
    env.update(OMP_NUM_THREADS=threads, MKL_NUM_THREADS=threads, OPENBLAS_NUM_THREADS="1")
    print(f"{directory.parent.name}/{directory.name}  {label}", flush=True)
    with log.open("wb") as f:
        run = subprocess.run(
            cmd,
            cwd=directory,
            stdout=f,
            stderr=subprocess.STDOUT,
            env=env,
            stdin=subprocess.DEVNULL,
        )
    output = log.read_text(encoding="utf-8", errors="replace")
    failures = (
        "convergence NOT achieved",
        "convergence NOT obtained",
        "Error in routine",
        "No convergence has been achieved",
        "Maximum CPU time exceeded",
    )
    ok = (
        run.returncode == 0
        and "JOB DONE." in output
        and not any(s.lower() in output.lower() for s in failures)
    )
    if kind == "relax":
        ok = ok and "bfgs converged" in output.lower()
    if not ok or not all((directory / p).is_file() for p in needed):
        raise RuntimeError(f"{label} nije potvrđen (kod {run.returncode}). Pogledaj {log}.")
    record = {"input_sha256": fingerprint, "command": [Path(x).name if x == binary else x for x in cmd]}
    if actual_input != contents:
        # Potvrda pamti stvarni ulaz i prepoznaje isti fizikalni račun u oba načina pohrane.
        record["workflow_sha256"] = fingerprint
        record["input_sha256"] = hashlib.sha256((actual_input + binary_signature).encode()).hexdigest()
    save_json(receipt, record)
    return log


def final_geometry(directory, g):
    tree = ET.parse(directory / "tmp/tise2.save/data-file-schema.xml")
    # QE zapisuje završne kartezijske koordinate u Bohrima. Zadržavamo isti redoslijed atoma.
    out = next(el for el in tree.getroot() if el.tag.split("}")[-1] == "output")
    block = out.find("atomic_structure/atomic_positions")
    atoms = list(block)
    if [el.get("name") for el in atoms] != g["symbols"]:
        raise ValueError("Redoslijed atoma nakon relaksacije ne odgovara referenci.")
    pos = np.array([[float(x) for x in el.text.split()] for el in atoms]) * BOHR_A
    return geom(g["cell_A"], g["symbols"], pos)


def energy(log):
    vals = re.findall(
        r"!\s+total energy\s*=\s*(" + REAL + r")\s+Ry", log.read_text(errors="replace")
    )
    if not vals:
        raise ValueError(f"QE energija nedostaje u {log}")
    return float(vals[-1].replace("D", "E")) * RY_EV

"""Dijagonalizacija, savijanje zone i praćenje modova. Frekvencije su u THz."""

import csv
import io
import re

import numpy as np
from scipy.constants import c as LIGHT_SPEED
from scipy.optimize import linear_sum_assignment

from qe import COPIES, MASSES, REAL, atomic_text, path_points


def csv_out(path, header, rows):
    buf = io.StringIO(newline="")
    out = csv.writer(buf, lineterminator="\n")
    out.writerow(header)
    out.writerows(rows)
    atomic_text(path, buf.getvalue())


def numbers(line):
    return [float(s.replace("D", "E").replace("d", "e")) for s in re.findall(REAL, line)]


def eigenfile(path, nat):
    lines = path.read_text(errors="replace").splitlines()
    freqs, modes = [], []
    for i, line in enumerate(lines):
        m = re.search(r"(?:freq|omega)\s*\(\s*\d+\s*\)\s*=\s*(" + REAL + r")\s*\[THz\]", line, re.I)
        if m:
            pairs = np.array([numbers(s) for s in lines[i + 1 : i + 1 + nat]])
            if pairs.shape != (nat, 6):
                raise ValueError(f"Neispravan oblik vektora u {path}.")
            modes.append((pairs[:, ::2] + 1j * pairs[:, 1::2]).ravel())
            freqs.append(float(m[1]))
    e = np.array(modes).T
    if e.shape != (3 * nat, 3 * nat) or not np.allclose(e.conj().T @ e, np.eye(3 * nat), atol=4e-5):
        raise ValueError(f"Očekuju se ortonormirani vektori iz fileig u {path}.")
    return np.array(freqs), e


def signed_sqrt(lam):
    return np.sign(lam) * np.sqrt(np.abs(lam))


def real_dynamical(nu, e):
    d = (e * (nu * np.abs(nu))) @ e.conj().T
    if np.max(np.abs(d.imag)) > 1e-4 * max(1.0, np.linalg.norm(d)):
        raise ValueError("D(Γ) ima neočekivano velik imaginarni dio.")
    return (d.real + d.real.T) / 2


def masses(g):
    return np.array([MASSES[s] for s in g["symbols"]])


def cartesian_direction(e, g):
    u = e.reshape(-1, 3) / np.sqrt(masses(g))[:, None]
    return u / np.linalg.norm(u, axis=1).max()


def orient(vec):
    return vec if vec[np.argmax(np.abs(vec))] >= 0 else -vec


SECTORS = {"Gamma": (0, 0), "M1": (1, 0), "M2": (0, 1), "M3": (1, 1)}
M_ARMS = ("M1", "M2", "M3")


def sector_modes(d):
    # Gamma superćelije sadrži četiri 9-dimenzionalna potprostora.
    modes = {}
    for name, (h, k) in SECTORS.items():
        signs = np.array([(-1.0) ** (h * i + k * j) for i, j in COPIES])
        b = np.kron(signs[:, None] / 2, np.eye(9))
        lam, e = np.linalg.eigh(b.T @ d @ b)
        modes[name] = signed_sqrt(lam), b @ e
    return modes


def frozen_directions(d, g, threshold, previous=None):
    """Izdvoji grane po translacijama i izgradi isti jednako težinski zbroj kao u bilježnici."""
    sectors = sector_modes(d)
    # Ako je na prvoj T samo neki drugi sektor nestabilan, sačuvaj i njegov zaseban PES.
    extra = previous["dodatni_sektor"] if previous else None
    if previous is None and sectors["M1"][0][0] >= -threshold:
        lowest = min(sectors, key=lambda label: sectors[label][0][0])
        if sectors[lowest][0][0] < -threshold:
            extra = lowest

    def describe(e):
        de = d @ e
        lam = float(e @ de)  # D je spremljena u predznačenim THz².
        return {
            "lambda_THz2": lam,
            "rezidual_THz2": float(np.linalg.norm(de - lam * e)),
            "negativna_zakrivljenost_iznad_praga": lam < -(threshold**2),
            "e_massweighted": e.tolist(),
            "u_rel": cartesian_direction(e, g).tolist(),
        }

    directions = {}
    names = list(M_ARMS)
    if extra and extra not in names:
        names.append(extra)
    for name in names:
        freq, modes = sectors[name]
        index, overlap = 0, None
        if previous is None:
            vec = modes[:, index]
            # Jedna primitivna ćelija ima devet komponenti; ostale kopije nose fazne znakove.
            vec = vec if vec[np.argmax(np.abs(vec[:9]))] >= 0 else -vec
        else:
            dots = modes.T @ previous["smjerovi"][name]["e_massweighted"]
            index = int(np.argmax(np.abs(dots)))
            overlap = float(abs(dots[index]))
            vec = modes[:, index] * (1 if dots[index] >= 0 else -1)
        directions[name] = {
            **describe(vec),
            "grana_sektora": index + 1,
            "nu_THz": float(freq[index]),
            "nestabilan": bool(freq[index] < -threshold),
            "preklop_s_prethodnom_T": overlap,
        }
    arms = np.column_stack([directions[name]["e_massweighted"] for name in M_ARMS])
    if not np.allclose(arms.T @ arms, np.eye(3), atol=1e-8):
        raise ValueError("Odabrana tri M-vektora moraju biti ortonormirana.")
    directions["3Q"] = {
        **describe(arms.sum(axis=1) / np.sqrt(3)),
        "koeficijenti_M1_M2_M3": [1 / np.sqrt(3)] * 3,
        "raspon_lambda_M_THz2": float(np.ptp([directions[k]["lambda_THz2"] for k in M_ARMS])),
    }
    return {"smjerovi": directions, "dodatni_sektor": extra}


def super_gamma(directory, g, threshold, reference=False):
    freq, raw = eigenfile(directory / "gamma.eig", len(g["symbols"]))
    d = real_dynamical(freq, raw)
    lam, modes = np.linalg.eigh(d)
    nu = signed_sqrt(lam)
    # Vanjski zapis čuva i izvorne QE vektore s proizvoljnom kompleksnom fazom.
    np.savez_compressed(
        directory / "gamma.npz",
        nu_THz=nu,
        e_massweighted=modes,
        qe_nu_THz=freq,
        qe_e_complex=raw,
        D_signed_THz2=d,
    )
    csv_out(
        directory / "gamma.csv",
        ["mod", "nu_THz", "nu2_THz2", "nestabilan"],
        [(i + 1, f, f * abs(f), bool(f < -threshold)) for i, f in enumerate(nu)],
    )
    rows = []
    for i in np.flatnonzero(nu < -threshold):
        v = orient(modes[:, i])
        u = cartesian_direction(v, g)
        for j, (s, vec, shift) in enumerate(zip(g["symbols"], v.reshape(-1, 3), u)):
            rows.append((i + 1, nu[i], j + 1, s, *vec, *shift))
    csv_out(
        directory / "nestabilni_vektori.csv",
        ["mod", "nu_THz", "atom", "element", "e_x", "e_y", "e_z", "u_x_rel", "u_y_rel", "u_z_rel"],
        rows,
    )
    if reference:
        sectors = sector_modes(d)
        reps = np.column_stack([ev for _, ev in sectors.values()])
        expected = np.concatenate([f * np.abs(f) for f, _ in sectors.values()])
        err = np.linalg.norm(d @ reps - reps * expected) / max(1.0, np.linalg.norm(d))
        if err > 0.005:
            raise ValueError(f"2×2 referenca nema očekivanu translacijsku simetriju ({err:.3g}).")
        csv_out(
            directory / "sektori.csv",
            ["sektor", "grana", "nu_THz", "relativna_pogreska_blokova"],
            [(s, i + 1, f, err) for s, (fs, _) in sectors.items() for i, f in enumerate(fs)],
        )
    return d


def frequency_file(path):
    text = path.read_text()
    header, data = text.split("/", 1)
    nb = int(re.search(r"nbnd\s*=\s*(\d+)", header)[1])
    nq = int(re.search(r"nks\s*=\s*(\d+)", header)[1])
    rows = iter(line for line in data.splitlines() if line.strip())
    points, frequencies = [], []
    for _ in range(nq):
        coord = numbers(next(rows))
        if len(coord) not in (3, 4):
            raise ValueError(f"Neočekivan redak q-točke u {path}.")
        points.append(coord[:3])  # QE 7.5 može dodati težinu kao četvrti broj.
        f = []
        while len(f) < nb:
            f.extend(numbers(next(rows)))
        if len(f) != nb:
            raise ValueError(f"Pogrešan broj frekvencija u {path}.")
        frequencies.append(f)
    if next(rows, None) is not None:
        raise ValueError(f"Višak redaka u {path}.")
    return np.array(points), np.array(frequencies) * (LIGHT_SPEED * 100 / 1e12)


def dispersion(directory, c, g):
    q = path_points(c)
    cart, nu = frequency_file(directory / "path.freq")
    if nu.shape != (len(q), 9):
        raise ValueError("Broj q-točaka ili grana ne odgovara primitivnoj ćeliji.")
    recip = 2 * np.pi * np.linalg.inv(g["cell_A"]).T
    alat = np.linalg.norm(g["cell_A"][0])
    if not np.allclose(cart, q @ recip * alat / (2 * np.pi), atol=2e-6):
        raise ValueError("Kartezijske q-točke iz QE ne odgovaraju zadanoj putanji.")
    dist = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(q @ recip, axis=0), axis=1))]
    csv_out(
        directory / "disperzija.csv",
        ["s_A^-1", "q1", "q2", "q3"] + [f"nu_{i}_THz" for i in range(1, 10)],
        [(s, *qi, *fi) for s, qi, fi in zip(dist, q, nu)],
    )
    ix, band = np.unravel_index(nu.argmin(), nu.shape)
    m = c["path_steps_per_segment"]
    unstable = np.argwhere(nu < -c["instability_threshold_THz"])
    csv_out(
        directory / "nestabilne_tocke.csv",
        ["q1", "q2", "q3", "grana", "nu_THz"],
        [(*q[i], j + 1, nu[i, j]) for i, j in unstable],
    )
    return {
        "q_min": q[ix].tolist(),
        "grana_min": int(band + 1),
        "nu_min_THz": float(nu[ix, band]),
        "nu_M_THz": float(nu[m].min()),
        "broj_nestabilnih_parova_q_grana": len(unstable),
    }


def optical(nu, modes, g):
    trans = np.kron(np.sqrt(masses(g))[:, None], np.eye(3))
    trans /= np.linalg.norm(trans, axis=0)
    # Projekcija na cijeli translacijski potprostor ne ovisi o izboru osi.
    score = np.sum((modes.T @ trans) ** 2, axis=1)
    indices = np.argsort(score)[-3:]
    if score[indices].sum() < 2.95:
        raise ValueError("Provjeri akustično pravilo sume i translacije u Γ.")
    keep = np.ones(len(nu), bool)
    keep[indices] = False
    return nu[keep], modes[:, keep]


def match(previous, cur_freq, cur_modes):
    overlaps = np.abs(previous.T @ cur_modes)
    rows, perm = linear_sum_assignment(-(overlaps**2))
    assert np.array_equal(rows, np.arange(len(perm)))
    modes = cur_modes[:, perm]
    modes *= np.where(np.sum(previous * modes, axis=0) < 0, -1, 1)
    return cur_freq[perm], modes, overlaps[rows, perm]

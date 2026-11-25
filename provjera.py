"""Brze provjere fizike i formata. Podaci su analitički, bez upisa u rezultate TiSe2."""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from zipfile import ZipFile

import numpy as np
from scipy.linalg import null_space

import modovi as m
import qe
import studija
from studija import settings


def main():
    c = settings(qe.ROOT / "postavke.json")
    qe.verify_pseudo()
    g = qe.primitive(c)
    sup = qe.supercell(g)
    dis = qe.distorted(sup, c)
    shift = np.array(dis["positions_A"]) - np.array(sup["positions_A"])
    assert np.allclose(np.average(shift, weights=m.masses(sup), axis=0), 0)
    assert np.isclose(np.sqrt(np.sum(shift**2) / 12), c["random_rms_A"])
    assert np.isclose(qe.KB_EV * 100 / qe.RY_EV, 0.000633362, atol=1e-9)
    trans = np.kron(np.sqrt(m.masses(g))[:, None], np.eye(3))
    trans /= np.linalg.norm(trans, axis=0)
    local = np.column_stack([trans, null_space(trans.T)])
    d = np.zeros((36, 36))
    for key, (h, k) in m.SECTORS.items():
        signs = np.array([(-1.0) ** (h * i + k * j) for i, j in qe.COPIES])
        basis = np.kron(signs[:, None] / 2, local)
        lam = np.arange(1.0, 10.0)
        if key == "Gamma":
            lam[:3] = 0.0
        else:
            lam[0] = -4.0
        d += (basis * lam) @ basis.T
    lam, e = np.linalg.eigh(d)
    nu = m.signed_sqrt(lam)
    # Slučajna kompleksna kombinacija triju degeneriranih M vektora.
    rng = np.random.default_rng(31)
    rotation = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))[0]
    raw = e.astype(complex)
    raw[:, :3] = raw[:, :3] @ rotation
    assert np.allclose(m.real_dynamical(nu, raw), d)
    sectors = m.sector_modes(d)
    for key in ("M1", "M2", "M3"):
        assert np.isclose(sectors[key][0][0], -2)
    u = m.cartesian_direction(sectors["M1"][1][:, 0], sup).reshape(4, 3, 3)
    assert np.allclose(u[0], u[1]) and np.allclose(u[0], -u[2])
    f, opt = m.optical(nu, e, sup)
    perm = rng.permutation(len(f))
    recovered, modes, overlap = m.match(opt, f[perm], -opt[:, perm])
    assert len(f) == 33 and np.allclose(recovered, f) and np.allclose(overlap, 1)
    with TemporaryDirectory() as tmp:
        p = Path(tmp)
        rows = []
        for i, fi in enumerate(nu):
            rows.append(f"freq ({i + 1:5}) = {fi:15.6f} [THz]")
            for vec in raw[:, i].reshape(12, 3):
                rows.append("( " + " ".join(f"{x.real:.6f} {x.imag:.6f}" for x in vec) + " )")
        (p / "gamma.eig").write_text("\n".join(rows))
        measured = m.super_gamma(p, sup, c["instability_threshold_THz"], reference=True)
        assert np.allclose(measured, d, atol=4e-5)
        assert len((p / "nestabilni_vektori.csv").read_text().splitlines()) == 1 + 3 * 12
        with patch.object(qe, "run_root", return_value=p):
            frozen_test(c, g, sup, d)
            tracking_test(c, g, sup, e, nu, p)
        pts = qe.path_points(c)
        cart = pts @ np.linalg.inv(g["cell_A"]).T * c["a_A"]
        for weight in ("", " 0.000000"):
            rows = [f" &plot nbnd=9, nks={len(pts)} /"]
            for q in cart:
                rows.extend([" ".join(str(x) for x in q) + weight, "-20 0 3 4 5 6", "7 8 9"])
            (p / "path.freq").write_text("\n".join(rows))
            result = m.dispersion(p, c, g)
            assert np.isclose(result["nu_M_THz"], -20 * 0.0299792458)
        lesson = qe.ROOT.parent / "vjezbe6.zip"
        if lesson.exists():
            with ZipFile(lesson) as z:
                (p / "si.modes").write_bytes(z.read("01_Si_Gamma/matdyn.modes"))
            # Jednaki atomi Si znače da su ovdje normirani pomaci također ortogonalni.
            f, _ = m.eigenfile(p / "si.modes", 2)
            assert len(f) == 6 and np.max(f) > 10
    print("OK: temperature, geometrija, tri M-sektora, translacije, preklopi i QE formati.")


def tracking_test(c, g, sup, e, nu, tmp):
    """Grana promijeni poredak, ostali se vektori na drugoj temperaturi permutiraju."""
    from izvjestaj import rel_branches

    test_config = dict(c, temperatures_K=[100, 1600])
    shuffle = np.random.default_rng(73).permutation(len(nu))
    for t in test_config["temperatures_K"]:
        p = studija.tdir(test_config, t)
        qe.save_json(p / "relax_1x1" / "geometrija.json", g)
        qe.save_json(p / "relax_2x2" / "geometrija.json", qe.distorted(sup, c))
        d = p / "gamma_2x2_rel"
        qe.save_json(d / "zavrseno.json", {})
        f = nu.copy()
        if t == 1600:
            f[0] = 4.0
            np.savez(d / "gamma.npz", nu_THz=f[shuffle], e_massweighted=e[:, shuffle])
        else:
            np.savez(d / "gamma.npz", nu_THz=f, e_massweighted=e)
    dest = tmp / "test_izvjestaja"
    dest.mkdir()
    rel_branches(test_config, dest)
    result = studija.read_json(dest / "najveci_pomak.json")
    assert np.isclose(result["raspon_THz"], 6) and result["sve_temperature"]
    assert np.isclose(result["najmanji_preklop"], 1)


def frozen_test(c, g, sup, d):
    """Dva različita presjeka poznatog Hessiana, ista referenca i proizvoljne faze iz eigh."""
    from scipy.constants import atomic_mass, elementary_charge

    t = c["temperatures_K"][0]
    sectors = m.sector_modes(d)
    v2, v3 = (sectors[key][1][:, 0] for key in ("M2", "M3"))
    split = d + 0.7 * np.outer(v2, v2) - 0.2 * np.outer(v3, v3)
    info = m.frozen_directions(split, sup, c["instability_threshold_THz"])
    arms = np.array([info["smjerovi"][name]["e_massweighted"] for name in m.M_ARMS])
    triple = info["smjerovi"]["3Q"]
    assert np.allclose(triple["e_massweighted"], arms.sum(axis=0) / np.sqrt(3))
    assert np.isclose(triple["lambda_THz2"], (-4 - 3.3 - 4.2) / 3)
    assert triple["rezidual_THz2"] > 0.1
    equal = m.frozen_directions(d, sup, c["instability_threshold_THz"])
    assert equal["smjerovi"]["3Q"]["rezidual_THz2"] < 1e-10
    crossing = split + 6.5 * np.outer(arms[0], arms[0])
    tracked = m.frozen_directions(crossing, sup, c["instability_threshold_THz"], info)
    m1 = tracked["smjerovi"]["M1"]
    assert m1["grana_sektora"] == 2 and m1["nu_THz"] > 0
    assert np.allclose(m1["e_massweighted"], arms[0])
    base = studija.tdir(c, t)
    qe.save_json(base / "gamma_2x2_ref" / "zavrseno.json", {"F_eV": -10})
    factor = atomic_mass * (2 * np.pi * 1e12) ** 2 * 1e-20 / elementary_charge

    def harmonic(directory, label, program, input_text):
        coords = input_text.split("ATOMIC_POSITIONS crystal\n")[1].split("K_POINTS")[0].splitlines()
        pos = np.array([[float(v) for v in line.split()[1:4]] for line in coords]) @ sup["cell_A"]
        weighted = ((pos - sup["positions_A"]) * np.sqrt(m.masses(sup))[:, None]).ravel()
        energy = -10 + 0.5 * factor * (weighted @ split @ weighted)
        log = directory / (label + ".out")
        qe.atomic_text(log, f"!    total energy = {energy / qe.RY_EV:.14f} Ry\n")
        return log

    with (
        patch.object(studija, "reference", return_value=split),
        patch.object(qe, "step", side_effect=harmonic) as scf,
    ):
        studija.pes(c, t, g)
        assert scf.call_count == 2 * (len(c["amplitudes_A"]) - 1)
    for name in ("M1", "3Q"):
        curve = base / "pes" / name
        stats = studija.read_json(curve / "zakrivljenost.json")
        data = np.loadtxt(curve / "potencijal.csv", delimiter=",", skiprows=1)
        assert stats["smjer"] == name and stats["minimum_na_rubu_intervala"]
        assert np.allclose(data[data[:, 0] == 0, 2:5], [[-10, 0, 0]])
        assert np.allclose(data[:, 3], data[:, 5], atol=1e-7)
        assert np.isclose(stats["PES_eV_A2"], stats["DFPT_eV_A2"], atol=1e-7)


if __name__ == "__main__":
    main()

"""Osjetljivost najnižeg moda M na k-mrežu, odsječenje ravnih valova i vakuum."""

import argparse
from pathlib import Path

import qe
from modovi import csv_out
from studija import prepare, read_json, run, settings, tdir


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", action="store_true", help="Izvrši DFPT za odabrane provjere.")
    p.add_argument("--temperatures", nargs="+", type=int, default=[100, 1600])
    p.add_argument("--config", type=Path, default=qe.ROOT / "postavke.json")
    a = p.parse_args()
    base = settings(a.config)
    if not set(a.temperatures) <= set(base["temperatures_K"]):
        p.error("Temperatura mora postojati u postavkama.")
    variants = {
        "osnovno": {},
        "k36": {"k_mesh": 36},
        "k48": {"k_mesh": 48},
        "odsijecanje": {"ecutwfc_Ry": 75, "ecutrho_Ry": 750},
        "vakuum": {"c_A": base["c_A"] + 6},
    }
    root = qe.run_root(base)
    rows = []
    for name, overrides in variants.items():
        config = dict(base, **overrides)
        prepare(config)
        for t in a.temperatures:
            if a.run:
                qe.verify_pseudo()
                run(config, t, "m", 1)
            f = tdir(config, t) / "fononi_1x1" / "M_direktno.json"
            value = read_json(f)["nu_min_THz"] if f.exists() else None
            rows.append(
                (
                    name,
                    t,
                    config["k_mesh"],
                    config["ecutwfc_Ry"],
                    config["ecutrho_Ry"],
                    config["c_A"],
                    value,
                )
            )
    csv_out(
        root / "rezultati" / "konvergencija_M.csv",
        ["varijanta", "T_K", "k", "ecutwfc_Ry", "ecutrho_Ry", "c_A", "nu_M_THz"],
        rows,
    )
    print(root / "rezultati" / "konvergencija_M.csv")


if __name__ == "__main__":
    main()

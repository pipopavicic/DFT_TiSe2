"""Preuzmi provjereno prijenosno izdanje QE 7.5 samo u mapu engine/."""

import hashlib
import os
from pathlib import Path
import shutil
import tempfile
from urllib.request import urlopen
from zipfile import ZipFile

BASE = "https://github.com/QMatSuite/quantum-espresso-windows-exe/releases/download/qe-7.5-win-oneapi-msmpi/"
ARCHIVE = "qe-7.5-win-oneapi-msmpi.zip"
SHA256 = "e50ec7368a5b7a964eb5084e0a13464ba2628d4e147c16227e876002ea34fdf7"
NAMES = {"pw.exe", "ph.exe", "q2r.exe", "matdyn.exe", "dynmat.exe", "mpiexec.exe"}


def main():
    if os.name != "nt":
        raise SystemExit("Za Linux postavi QE_BIN ili dodaj programe Quantum ESPRESSO u PATH.")
    dest = Path(__file__).resolve().parent / "engine"
    dest.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        file = Path(tmp) / ARCHIVE
        print("Preuzimam QE 7.5 (400 MB).", flush=True)
        with urlopen(BASE + ARCHIVE, timeout=90) as src, file.open("wb") as out:
            shutil.copyfileobj(src, out, 1024 * 1024)
        h = hashlib.sha256()
        with file.open("rb") as src:
            for block in iter(lambda: src.read(1024 * 1024), b""):
                h.update(block)
        if h.hexdigest() != SHA256:
            raise ValueError("SHA-256 arhive ne odgovara objavljenom izdanju.")
        with ZipFile(file) as z:
            for item in z.namelist():
                name = Path(item).name
                if name in NAMES or name.lower().endswith(".dll"):
                    with z.open(item) as src, (dest / name).open("wb") as out:
                        shutil.copyfileobj(src, out)
    print(dest)


if __name__ == "__main__":
    main()

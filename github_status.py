"""Na Srcu prati već pokrenutu DFT studiju i ažuriraj njezin GitHub Issue."""

import argparse
from datetime import datetime, timezone
import getpass
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import warnings
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError


ISSUE = "https://github.com/pipopavicic/DFT_TiSe2/issues/1"
API = "https://api.github.com/repos/pipopavicic/DFT_TiSe2/issues/1"
MARKER = "<!-- tise2-dft-status -->"
INTERVAL = 20 * 60
IDLE_GRACE = 90
POLL = 15
STAGES = (
    # Iste izlazne datoteke koristi i `studija.py status` preko modula izvjestaj.
    ("Relaks. 1×1", "relax_1x1", "geometrija.json", "relax.out"),
    ("M", "fononi_1x1", "M_direktno.json", "M.out"),
    ("Disperzija", "fononi_1x1", "minimum.json", "mesh.out"),
    ("Γ ref.", "gamma_2x2_ref", "zavrseno.json", "scf.out"),
    ("PES M1", "pes/M1", "zakrivljenost.json", "pomak_*/scf.out"),
    ("PES 3Q", "pes/3Q", "zakrivljenost.json", "pomak_*/scf.out"),
    ("Relaks. 2×2", "relax_2x2", "geometrija.json", "relax.out"),
    ("Γ rel.", "gamma_2x2_rel", "zavrseno.json", "scf.out"),
)
QE = {"ph.x", "pw.x", "q2r.x", "matdyn.x", "dynmat.x"}
Q_TOTAL = re.compile(r"^\s*\(\s*(\d+)\s+q-points?\)")
LAST_MODE = re.compile(r"^\s*freq\s*\(\s*9\s*\)\s*=")


def q_progress(log):
    """Za svaki q primitivne ćelije TiSe2 brojimo zadnju od devet grana."""
    total, done = None, 0
    counted = False
    try:
        with log.open(encoding="utf-8", errors="replace") as src:
            for line in src:
                if total is None:
                    match = Q_TOTAL.match(line)
                    if match:
                        total = int(match[1])
                if "Calculation of q =" in line:
                    counted = False
                if not counted and LAST_MODE.match(line):
                    done += 1
                    counted = True
    except FileNotFoundError:
        pass
    return f"{done}/{total} q" if total else ""


def stage_name(parts, input_name):
    folder = parts[0]
    if folder == "fononi_1x1":
        return {
            "scf.in": "SCF 1×1", "M.in": "Fonon u M", "mesh.in": "Disperzija 1×1",
            "q2r.in": "Konstante sila 1×1", "matdyn.in": "Interpolacija disperzije",
        }.get(input_name, "Fononi 1×1")
    if folder == "pes" and len(parts) >= 2 and parts[1] in ("M1", "M2", "M3", "3Q"):
        return "PES " + parts[1]
    return {
        "relax_1x1": "Relaksacija 1×1", "relax_2x2": "Relaksacija 2×2",
        "gamma_2x2_ref": "Γ referentne 2×2", "gamma_2x2_rel": "Γ relaksirane 2×2",
    }.get(folder, "QE račun")


def holds_study_lock(proc, locks):
    # Otvoreni deskriptori otkrivaju aktivnu studiju i između dva poziva QE-u.
    # Samo postojanje njezine .lock datoteke nije znak aktivnosti.
    try:
        for fd in (proc / "fd").iterdir():
            try:
                if os.readlink(fd) in locks:
                    return True
            except OSError:
                continue
    except OSError:
        pass
    return False


def active_work(root, temperatures):
    names = {f"T{t:04d}K": t for t in temperatures}
    locks = {str(root / name / ".izvodenje.lock") for name in names}
    locks.add(str(root / ".cijela_serija.lock"))
    processes = subprocess.check_output(
        ["ps", "-u", str(os.getuid()), "-o", "pid=,etimes=,comm="], text=True
    )
    driver_alive, jobs = False, {}
    for row in processes.splitlines():
        pid, elapsed, command = row.split(maxsplit=2)
        proc = Path("/proc") / pid
        if command.startswith(("python", "ipython")) and not driver_alive:
            driver_alive = holds_study_lock(proc, locks)
        if command not in QE:
            continue
        try:
            parts = Path(os.readlink(proc / "cwd")).relative_to(root).parts
            args = (proc / "cmdline").read_bytes().split(b"\0")
        except (OSError, ValueError):
            continue
        if len(parts) < 2 or parts[0] not in names:
            continue
        input_name = next((os.fsdecode(arg) for arg in args if arg.endswith(b".in")), "")
        key = (names[parts[0]], stage_name(parts[1:], Path(input_name).name))
        jobs[key] = max(jobs.get(key, 0), int(elapsed))
    return driver_alive or bool(jobs), jobs


def snapshot(root, temperatures, jobs, message):
    now = datetime.now(timezone.utc)
    try:
        now = now.astimezone(ZoneInfo("Europe/Zagreb"))
    except ZoneInfoNotFoundError:
        pass
    disk = shutil.disk_usage(root)
    lines = [MARKER, f"**{message}**", "", f"Očitano {now:%d.%m.%Y. %H:%M %Z}.",
             f"Slobodno na disku s rezultatima: {disk.free / 2**30:.1f} GiB.", ""]
    for (temp, phase), seconds in sorted(jobs.items()):
        lines.append(f"- **{temp} K** · {phase} · proces traje {seconds // 3600} h "
                     f"{seconds % 3600 // 60} min")
    if jobs:
        lines.append("")
    lines += ["| T (K) | " + " | ".join(stage[0] for stage in STAGES) + " |",
              "|" + "---|" * (len(STAGES) + 1)]
    for t in temperatures:
        temperature_dir = root / f"T{t:04d}K"
        cells = [str(t)]
        for _, folder, marker, output in STAGES:
            directory = temperature_dir / folder
            if (directory / marker).is_file():
                value = "gotovo"
            elif any(directory.glob(output)):
                value = "započeto"
            else:
                value = "čeka"
            if output == "mesh.out":
                progress = q_progress(directory / output)
                if progress:
                    value += f" ({progress})"
            cells.append(value)
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "Faze se očitavaju kao u `studija.py status`. Q-točka je dovršena kad "
              "Quantum ESPRESSO za nju ispiše svih devet frekvencija primitivne ćelije."]
    return "\n".join(lines)


def publish(body, token):
    req = Request(API, method="PATCH", data=json.dumps({"body": body}).encode(),
                  headers={"Accept": "application/vnd.github+json", "Content-Type": "application/json",
                           "User-Agent": "DFT-TiSe2-status"})
    req.add_unredirected_header("Authorization", "Bearer " + token)
    try:
        with urlopen(req, timeout=15):
            return True
    except HTTPError as exc:
        limited = exc.headers.get("Retry-After") or exc.headers.get("X-RateLimit-Remaining") == "0"
        if exc.code < 500 and exc.code != 429 and not limited:
            raise SystemExit(f"GitHub je odbio upis (HTTP {exc.code}). Provjeri rok tokena, "
                             "repozitorij DFT_TiSe2 i dozvolu Issues → Read and write.") from None
        print(f"GitHub privremeno nedostupan (HTTP {exc.code}).", flush=True)
    except (URLError, OSError):
        print("Mrežna veza prema GitHubu nije uspjela.", flush=True)
    return False


def monitor(root, temperatures, token):
    last_active, next_send = time.monotonic(), 0
    while True:
        active, jobs = active_work(root, temperatures)
        now = time.monotonic()
        if active:
            last_active = now
        finished = now - last_active >= IDLE_GRACE
        if now >= next_send or finished:
            if active:
                message = "DFT račun je aktivan."
            elif finished:
                message = "Praćenje završeno. Proces studije više nije aktivan."
            else:
                message = "Trenutačno nema aktivnog koraka. Čekam sljedeći proces studije."
            body = snapshot(root, temperatures, jobs, message)
            print(body.split("\n", 1)[1], flush=True)
            success = publish(body, token)
            if success:
                print("Objavljeno:", ISSUE, flush=True)
            if finished:
                if not success:
                    raise SystemExit("Završno očitanje je iznad. Objava na GitHubu nije uspjela.")
                return
            next_send = now + INTERVAL
        time.sleep(POLL)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("racun", type=Path, help="Npr. racuni/k24_q6_17503ad3bf")
    parser.add_argument("--preview", action="store_true", help="Jedno lokalno očitanje")
    args = parser.parse_args()
    root = args.racun.resolve()
    try:
        temperatures = json.loads((root / "postavke.json").read_text(encoding="utf-8-sig"))["temperatures_K"]
    except (OSError, KeyError, ValueError) as exc:
        parser.error(f"Ne mogu učitati temperature iz {root / 'postavke.json'} ({exc}).")
    if not isinstance(temperatures, list) or not temperatures or not all(
        type(t) is int and t > 0 for t in temperatures
    ):
        parser.error("Očekujem popis pozitivnih cjelobrojnih temperatura.")
    if args.preview:
        if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(encoding="utf-8")
        jobs = active_work(root, temperatures)[1] if sys.platform.startswith("linux") else {}
        print(snapshot(root, temperatures, jobs, "Lokalni pregled."))
        return
    if not sys.platform.startswith("linux"):
        parser.error("Prati procese pokretanjem skripte u Jupyteru na Srcu.")
    import fcntl

    with (root.parent / ".github_status.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("Jedno praćenje već radi u ovom projektu.") from None
        if not active_work(root, temperatures)[0]:
            raise SystemExit("Nema aktivnog računa u odabranoj mapi. Pokreni praćenje dok račun radi.")
        token = os.environ.get("TISE2_GITHUB_TOKEN")
        if not token:
            with warnings.catch_warnings():
                warnings.simplefilter("error", getpass.GetPassWarning)
                try:
                    token = getpass.getpass("GitHub token za status (skriveni unos): ").strip()
                except getpass.GetPassWarning:
                    raise SystemExit("Za skriveni unos pokreni skriptu preko %run u Python bilježnici.") from None
        if not token:
            raise SystemExit("Za objavu je potreban GitHub token.")
        monitor(root, temperatures, token)


if __name__ == "__main__":
    main()

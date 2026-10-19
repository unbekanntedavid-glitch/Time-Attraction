#!/usr/bin/env python3
"""
SPARC full runner – egy közös hibametrika, hat modell.
Relatív bemeneti/kimeneti útvonalak, bemeneti hash-ellenőrzés,
hiba hiányzó galaxis esetén.

Bemenet (relatív a szkript mappájához vagy aktuális könyvtárhoz):
  data/rotmod/*_rotmod.dat
  data/blind_123.txt

Kimenet:
  output/points_<model>.csv
  output/summary_<model>.csv

Rotmod-archívum beszerzése:
  A SPARC Rotmod_LTG.zip a http://astroweb.cwru.edu/SPARC/ oldalról
  vagy a Lelli et al. 2016 (AJ 152, 157) kiegészítő anyagából tölthető le.
  Csomagold ki a data/rotmod/ mappába úgy, hogy a fájlok
  data/rotmod/<galaxis>_rotmod.dat formában legyenek.
"""
import os, sys, csv, math, statistics, hashlib, zipfile, tempfile, shutil
from pathlib import Path

# --- relatív útvonalak ---
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
ROTMOD = DATA / "rotmod"
BLIND_FILE = DATA / "blind_123.txt"
OUTPUT = ROOT / "output"

SIGMA_FLOOR = 20.0
SIGMA_FRAC = 0.05

# Várt bemeneti hash (SHA256 a blind_123.txt-re; Rotmod fájlok egyenként ellenőrizhetők)
EXPECTED_BLIND_HASH = "3637d2f5bda6b2b97e645c5fba5cf39d310a1826bee388286c2c74b9171aa283"
EXPECTED_ROTMOD_HASH = "0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588"
EXPECTED_RADIUS_POINTS = 2596

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def sigma(Vobs):
    return max(SIGMA_FRAC * Vobs, SIGMA_FLOOR)

def load_blind():
    if not BLIND_FILE.exists():
        sys.exit(f"HIBA: hiányzó bemenet: {BLIND_FILE}")
    h = sha256_file(BLIND_FILE)
    print(f"blind_123.txt SHA256: {h}")
    if EXPECTED_BLIND_HASH and h != EXPECTED_BLIND_HASH:
        sys.exit(f"HIBA: blind_123.txt hash eltérés (várt {EXPECTED_BLIND_HASH})")
    return set(BLIND_FILE.read_text().split())

def parse_rotmod(fp):
    rows = []
    with open(fp) as f:
        for line in f:
            p = line.split()
            if len(p) < 6:
                continue
            try:
                rows.append((float(p[0]), float(p[1]), float(p[3]), float(p[4]), float(p[5])))
            except ValueError:
                continue
    if rows and rows[0][0] < 0.01:
        rows[0] = (0.16,) + rows[0][1:]
    return rows

def run_model(name, vpred_fn, blind, outdir):
    points = []
    summary = []
    missing = []
    for gal in sorted(blind):
        fp = ROTMOD / f"{gal}_rotmod.dat"
        if not fp.exists():
            missing.append(gal)
            continue
        raw = parse_rotmod(fp)
        if not raw:
            missing.append(gal)
            continue
        chi2 = 0.0
        n = 0
        for R, Vobs, Vg, Vd, Vb in raw:
            Vbar = math.sqrt(Vg**2 + Vd**2 + Vb**2)
            Vpred = vpred_fn(R, Vbar, Vg, Vd, Vb)
            if Vobs > 0:
                s = sigma(Vobs)
                n += 1
                chi2 += ((Vpred - Vobs) / s)**2
                points.append({
                    "model": name, "galaxy": gal, "R_kpc": R,
                    "Vobs": Vobs, "Vbar": Vbar, "Vpred": Vpred,
                    "sigma": s, "residual": (Vpred - Vobs) / s
                })
        if n > 1:
            summary.append({"model": name, "galaxy": gal, "n_radii": n,
                            "chi2_red": chi2 / (n - 1)})
    if missing:
        sys.exit(f"HIBA [{name}]: hiányzó galaxisok ({len(missing)}): {missing[:10]}")
    pf = outdir / f"points_{name}.csv"
    with open(pf, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["model","galaxy","R_kpc","Vobs","Vbar","Vpred","sigma","residual"])
        w.writeheader()
        w.writerows(points)
    sf = outdir / f"summary_{name}.csv"
    with open(sf, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["model","galaxy","n_radii","chi2_red"])
        w.writeheader()
        w.writerows(summary)
    med = statistics.median(r["chi2_red"] for r in summary) if summary else float("nan")
    print(f"{name}: {len(summary)} galaxis, medián χ²_red = {med:.4f}")
    return med


def verify_rotmod_archive():
    """Check Rotmod zip hash if present, and count radius points in blind sample."""
    zip_path = ROOT / "data" / "Rotmod_LTG.zip"
    if not zip_path.exists():
        sys.exit("HIBA: hiányzó Rotmod_LTG.zip a data/ mappában (exit 1)")
    h = sha256_file(zip_path)
    print(f"Rotmod_LTG.zip SHA256: {h}")
    if h != EXPECTED_ROTMOD_HASH:
        sys.exit(f"HIBA: Rotmod hash eltérés (várt {EXPECTED_ROTMOD_HASH})")
    # Extract to a fresh empty directory after verification
    extract_dir = ROOT / "data" / "rotmod_extracted"
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True)
    with zipfile.ZipFile(zip_path) as z:
        z.extractall(extract_dir)
    print(f"Rotmod kicsomagolva: {extract_dir}")
    # Point ROTMOD to the extracted location
    global ROTMOD
    ROTMOD = extract_dir
    # count points
    n = 0
    for gal in load_blind():
        fp = ROTMOD / f"{gal}_rotmod.dat"
        if fp.exists():
            for line in fp.read_text().splitlines():
                p = line.split()
                if len(p) >= 6:
                    try:
                        float(p[0]); float(p[1])
                        n += 1
                    except: pass
    print(f"Blind minta sugárpontok: {n}")
    if n != EXPECTED_RADIUS_POINTS:
        sys.exit(f"HIBA: sugárpont szám eltérés (várt {EXPECTED_RADIUS_POINTS}, kapott {n})")

def main():
    verify_rotmod_archive()
    blind = load_blind()
    OUTPUT.mkdir(exist_ok=True)
    print(f"Bemenet: {ROTMOD}")
    print(f"Kimenet: {OUTPUT}")

    run_model("vbar", lambda R, Vb, *a: Vb, blind, OUTPUT)
    run_model("vbar_095", lambda R, Vb, *a: Vb * 0.95, blind, OUTPUT)
    run_model("vbar_105", lambda R, Vb, *a: Vb * 1.05, blind, OUTPUT)
    run_model("no_bulge", lambda R, Vb, Vg, Vd, Vbul: math.sqrt(Vg**2 + Vd**2), blind, OUTPUT)
    run_model("no_gas", lambda R, Vb, Vg, Vd, Vbul: math.sqrt(Vd**2 + Vbul**2), blind, OUTPUT)
    run_model("disk_only", lambda R, Vb, Vg, Vd, Vbul: abs(Vd), blind, OUTPUT)
    print("Kész.")

if __name__ == "__main__":
    main()

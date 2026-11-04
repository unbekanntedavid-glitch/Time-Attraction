# SPDX-License-Identifier: CC-BY-NC-4.0
# Copyright (c) 2026 Szabolcsi David (ORCID 0009-0005-6300-074X). Commercial use only with the author's written permission.
"""Check every numeric claim of the article against the module archives.

usage (from the article folder): python3 src/check_claims.py [modules_dir]     (default: .. , the set root)

  1  every claim is re-read from the archives; value, EN and HU strings equal CLAIMS_SOURCES.csv
  2  values_*.tex, tab_*.tex and the claim table are reproduced byte-for-byte by make_values.py
  3  every value macro used in the article exists; EN and HU use the same set; every claim is used
  4  no hand-typed number in the article body outside the identifier/structure whitelist (src/texscan.py)
  5  cross-module equalities (the same quantity read from two modules)
  6  the clock module's bound files are byte-identical to the upstream module outputs and match PROVENANCE_bounds.json
  7  the article figures are byte-identical to the module files
  8  only the project DOIs; every file name cited with \\file{} exists in the set or in a module archive
  9  the PDFs exist and carry the licence footer on every page
 10  the fresh-run logs contain the results quoted in the article
Exit code 0 only if every check passes.
"""
import csv
import hashlib
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from claimlib import MODULES, Store, evaluate, fmt, hu, plain  # noqa: E402
from claims_spec import CLAIMS, EQUAL  # noqa: E402
import make_values  # noqa: E402
import texscan  # noqa: E402
from build_article import FIGURES, PDF  # noqa: E402

DOIS = {"10.5281/zenodo.22983681", "10.5281/zenodo.23540741", "10.5281/zenodo.23727585", "10.5281/zenodo.23814053"}
BOUNDS = {  # clock inputs/bounds/<name> -> (upstream module, file)
    "bren_continuation_results_r2.json": ("bren", "outputs/bren_continuation_results.json"),
    "subcm_envelopes.npz": ("subcm", "outputs/envelopes.npz"),
    "crossings.csv": ("subcm", "outputs/crossings.csv"),
    "AHN2003_fig4__fill_polygons.csv": ("subcm", "outputs/digitized/AHN2003_fig4__fill_polygons.csv"),
}
results = []


def check(name, ok, detail=""):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))


def main(modules_dir):
    art = os.getcwd()
    st = Store(modules_dir, art)
    rows = {r["id"]: r for r in csv.DictReader(open("CLAIMS_SOURCES.csv", encoding="utf-8"))}

    # 1
    bad = []
    for c in CLAIMS:
        x = evaluate(c, st)
        r = rows.get(c["id"])
        if r is None:
            bad.append(f"{c['id']}: not in table"); continue
        if isinstance(x, str):
            same = r["value_raw"] == repr(x)
        else:
            y = float(r["value_raw"])
            same = math.isfinite(x) and (x == y or abs(x - y) <= 1e-12 * max(abs(x), abs(y)))
        en = fmt(x, c["fmt"]); hu_ = hu(en) if c["fmt"][0] != "str" else en
        if not same or plain(en) != r["printed_en"] or plain(hu_) != r["printed_hu"]:
            bad.append(c["id"])
    check(f"1 claims re-read from the archives ({len(CLAIMS)})", not bad and len(rows) == len(CLAIMS), ", ".join(bad[:10]))

    # 2
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copytree("tex", os.path.join(tmp, "tex"))
        shutil.copytree("logs", os.path.join(tmp, "logs"))
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            make_values.main(modules_dir, tmp)
        diff = []
        for rel in ["CLAIMS_SOURCES.csv", "CLAIMS_SOURCES.md"] + ["generated/" + f for f in sorted(os.listdir(os.path.join(tmp, "generated")))]:
            a = os.path.join(tmp, rel)
            if not os.path.exists(rel) or open(a, "rb").read() != open(rel, "rb").read():
                diff.append(rel)
        extra = set(os.listdir("generated")) - set(os.listdir(os.path.join(tmp, "generated")))
        diff += sorted("generated/" + e for e in extra)
    check("2 generated values, tables and claim table reproduced byte-for-byte", not diff, ", ".join(diff))

    # 3
    ids = {c["id"] for c in CLAIMS}
    use = {}
    for lang in ("en", "hu"):
        t = texscan.expand_tables(open(f"tex/article_{lang}.tex", encoding="utf-8").read(), "generated")
        use[lang] = set(texscan.macros(t))
    missing = (use["en"] | use["hu"]) - ids
    check("3a every value macro exists", not missing, ", ".join(sorted(missing)))
    check("3b EN and HU use the same value macros", use["en"] == use["hu"], ", ".join(sorted(use["en"] ^ use["hu"])))
    unused = ids - use["en"] - {a for p in EQUAL for a in p}
    check("3c every claim of the table is used in the article", not unused, ", ".join(sorted(unused)))
    sec_missing = [i for i in use["en"] if not rows[i]["sections_en"] or not rows[i]["sections_hu"]]
    check("3d every used claim has its sections recorded (EN, HU)", not sec_missing, ", ".join(sec_missing[:10]))

    # 4
    for lang in ("en", "hu"):
        t = texscan.expand_tables(texscan.body(f"tex/article_{lang}.tex"), "generated")
        h = texscan.hand_typed(t)
        check(f"4 no hand-typed numbers in the article body ({lang.upper()})", not h,
              "; ".join(f"line {i}: {tok!r}" for i, tok, _ in h[:8]))

    # 5
    bad = []
    for a, b in EQUAL:
        xa = float(rows[a]["value_raw"]); xb = float(rows[b]["value_raw"])
        if abs(xa - xb) > 1e-12 * max(abs(xa), abs(xb)):
            bad.append(f"{a}={xa} vs {b}={xb}")
    check(f"5 cross-module equalities ({len(EQUAL)})", not bad, "; ".join(bad))

    # 6
    import json
    prov = json.loads(st.raw("clock", "inputs/PROVENANCE_bounds.json"))["files"]
    bad = []
    for name, (mod, rel) in BOUNDS.items():
        local = st.raw("clock", "inputs/bounds/" + name)
        up = st.raw(mod, rel)
        sha = hashlib.sha256(local).hexdigest()
        if local != up:
            bad.append(f"{name} differs from {mod}:{rel}")
        if prov.get("bounds/" + name, {}).get("sha256") != sha:
            bad.append(f"{name} SHA not in PROVENANCE")
    if set(prov) != {"bounds/" + n for n in BOUNDS}:
        bad.append("PROVENANCE lists other files")
    check("6 bound files byte-identical to upstream and to PROVENANCE_bounds.json (4)", not bad, "; ".join(bad))

    # 7
    bad = [name for mod, rel, name in FIGURES
           if open(os.path.join("figures", name), "rb").read() != st.raw(mod, rel)]
    check(f"7 figures byte-identical to the module files ({len(FIGURES)})", not bad, ", ".join(bad))

    # 8
    names = set()
    for mod, (zname, prefix) in MODULES.items():
        for n in zipfile.ZipFile(os.path.join(modules_dir, zname)).namelist():
            parts = n[len(prefix):].rstrip("/").split("/")
            names.update(parts)
    names.update(os.listdir(modules_dir)); names.update(os.listdir(".")); names.add(os.path.basename(art))
    bad = []
    for lang in ("en", "hu"):
        t = open(f"tex/article_{lang}.tex", encoding="utf-8").read()
        for d in re.findall(r"\\doi\{([^}]*)\}", t) + re.findall(r"zenodo\.\d+", re.sub(r"\\doi\{[^}]*\}", "", t)):
            if d not in DOIS:
                bad.append(f"DOI {d}")
        for f in re.findall(r"\\file\{([^}]*)\}", t):
            if f.rstrip("/").split("/")[-1] not in names:
                bad.append(f"file {f}")
        for f in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]*)\}", t):
            if not os.path.exists(os.path.join("tex", f)):
                bad.append(f"figure {f}")
    check("8 only project DOIs; every cited file exists", not bad, "; ".join(sorted(set(bad))))

    # 9
    bad = []
    for lang, pdf in PDF.items():
        if not os.path.exists(pdf):
            bad.append(f"{pdf} missing"); continue
        info = subprocess.run(["pdfinfo", pdf], capture_output=True, text=True).stdout
        n = int(re.search(r"Pages:\s+(\d+)", info).group(1))
        txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True).stdout
        pages = txt.split("\f")[:n]
        key = "CC BY-NC 4.0 International" if lang == "en" else "CC BY-NC 4.0 Nemzetközi"
        miss = [i + 1 for i, p in enumerate(pages) if key not in p or "0009-0005-6300-074X" not in p]
        if miss:
            bad.append(f"{pdf}: no footer on pages {miss}")
        print(f"       {pdf}: {n} pages")
    check("9 PDFs present, licence footer on every page", not bad, "; ".join(bad))

    # 10
    bad = []
    for mod, line in (("clock", "RESULT: PASS"), ("bren", "all steps passed"), ("subcm", "all steps passed (frozen-digitization mode)")):
        log = open(os.path.join("logs", f"fresh_run_{mod}.txt"), encoding="utf-8").read()
        if line not in log or not re.search(r"^exit(?: code)?:?\s*0\s*$", log, re.I | re.M):
            bad.append(mod)
    check("10 fresh-run logs show the quoted results and exit code 0", not bad, ", ".join(bad))

    print(f"{sum(results)}/{len(results)} checks passed")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else ".."))

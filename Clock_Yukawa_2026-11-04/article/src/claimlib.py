# SPDX-License-Identifier: CC-BY-NC-4.0
# Copyright (c) 2026 Szabolcsi David (ORCID 0009-0005-6300-074X). Commercial use only with the author's written permission.
"""Read values from the module archives (read-only, no unpacking) and format them for the article."""
import csv
import io
import json
import math
import os
import re
import zipfile
from fractions import Fraction

MODULES = {
    "clock": ("tal_branch_clock_accel_test_r4.zip", "tal_branch_clock_accel_test/"),
    "bren": ("bren_continuation_zenodo_r4.zip", "bren_continuation/"),
    "subcm": ("subcm_yukawa_window_zenodo_r4.zip", "subcm_yukawa_window/"),
}


class Store:
    """Access to package files. `modules_dir` holds the three r4 zips; `article_dir` holds the article files."""

    def __init__(self, modules_dir, article_dir):
        self.modules_dir = modules_dir
        self.article_dir = article_dir
        self._zips = {}
        self._cache = {}

    def raw(self, module, rel):
        key = (module, rel)
        if key in self._cache:
            return self._cache[key]
        if module == "article":
            with open(os.path.join(self.article_dir, rel), "rb") as f:
                data = f.read()
        else:
            zname, prefix = MODULES[module]
            if module not in self._zips:
                self._zips[module] = zipfile.ZipFile(os.path.join(self.modules_dir, zname))
            data = self._zips[module].read(prefix + rel)
        self._cache[key] = data
        return data

    def get(self, src):
        """src = (module, file, kind, selector)."""
        module, rel, kind, sel = src
        data = self.raw(module, rel)
        if kind == "json":
            x = json.loads(data.decode("utf-8"))
            for k in sel:
                x = x[k]
            return x
        if kind == "csv":
            # sel = (filter dict, column); comment lines starting with '#' are skipped
            flt, col = sel
            lines = [l for l in data.decode("utf-8").splitlines() if not l.startswith("#")]
            rows = list(csv.DictReader(io.StringIO("\n".join(lines))))
            hit = [r for r in rows if all(r[k] == v for k, v in flt.items())]
            if col == "#count":
                return len(hit)
            if len(hit) != 1:
                raise KeyError(f"{rel}: {flt} matched {len(hit)} rows")
            return float(hit[0][col])
        if kind == "re":
            # sel = (pattern, group, cast)
            pat, grp, cast = sel
            m = re.search(pat, data.decode("utf-8"), re.M)
            if not m:
                raise KeyError(f"{rel}: pattern {pat!r} not found")
            v = m.group(grp)
            return {"float": float, "int": int, "str": str}[cast](v)
        if kind == "jkey":
            # sel = (path to a dict, key): the key must exist; its name is the value
            path, key = sel
            x = json.loads(data.decode("utf-8"))
            for k in path:
                x = x[k]
            if key not in x:
                raise KeyError(f"{rel}: key {key!r} not found")
            return float(key)
        if kind == "re_count":
            return len(re.findall(sel, data.decode("utf-8"), re.M))
        raise ValueError(kind)


def sel_str(src):
    module, rel, kind, sel = src
    if kind == "json":
        s = "/".join(str(k) for k in sel)
    elif kind == "csv":
        flt, col = sel
        s = "row " + ",".join(f"{k}={v}" for k, v in flt.items()) + f" col {col}"
    elif kind == "re":
        s = f"regex {sel[0]!r}"
    elif kind == "jkey":
        s = "key " + "/".join(str(k) for k in sel[0]) + f"/{sel[1]}"
    else:
        s = f"count {sel!r}"
    return s


def _sig_digits(x, n):
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    dec = max(0, n - 1 - e)
    return f"{x:.{dec}f}"


def fmt(x, spec):
    """Return the EN string (math content, no $). HU is derived by hu()."""
    kind = spec[0]
    if kind == "str":
        return str(x)
    if kind == "int":
        if float(x) != int(round(float(x))):
            raise ValueError(f"not an integer: {x}")
        return f"{int(round(float(x)))}"
    if kind == "fix":
        return f"{x:.{spec[1]}f}"
    if kind == "sig":
        return _sig_digits(x, spec[1])
    if kind == "sci":
        if x == 0:
            return "0"
        e = int(math.floor(math.log10(abs(x))))
        m = x / 10 ** e
        if round(abs(m), spec[1] - 1) >= 10:
            e += 1
            m = x / 10 ** e
        return f"{m:.{spec[1] - 1}f}\\times10^{{{e}}}"
    if kind == "frac":
        f = Fraction(x).limit_denominator(12)
        if abs(float(f) - x) > 1e-12:
            raise ValueError(f"not a simple fraction: {x}")
        return f"{f.numerator}" if f.denominator == 1 else f"{f.numerator}/{f.denominator}"
    raise ValueError(spec)


def hu(s):
    """Decimal point -> decimal comma (math mode); minus signs unchanged."""
    return re.sub(r"(\d)\.(\d)", r"\1{,}\2", s)


def plain(s):
    """Plain-text normal form of a formatted value (for the PDF number comparison and the CSV)."""
    s = s.replace("{,}", ",").replace("\\times10^{", "e").replace("}", "")
    return s


def evaluate(claim, store):
    vals = [store.get(s) for s in claim["src"]]
    expr = claim.get("expr", "a")
    ns = {"__builtins__": {}, "math": math, "abs": abs, "min": min, "max": max, "len": len, "sum": sum,
          "round": round, "float": float}
    for i, v in enumerate(vals):
        ns["abcdefgh"[i]] = v
    return eval(expr, ns)

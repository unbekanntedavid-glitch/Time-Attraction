#!/usr/bin/env python3
"""Inline data, numerical core, chart helper, i18n and UI into ../index.html."""
import os
here = os.path.dirname(os.path.abspath(__file__))
rd = lambda n: open(os.path.join(here, n), encoding="utf-8").read()
html = rd("template.html")
for tag, fn in [("CSS", "style.css"), ("DATA", "data.js"), ("CORE", "core.js"), ("CHART", "chart.js"), ("I18N", "i18n.js"), ("REPRO", "repro.js"), ("UI", "ui.js")]:
    body = rd(fn)
    if tag != "CSS": body = body.replace("</script", "<\\/script")
    html = html.replace("/*@@%s@@*/" % tag, body)
out = os.path.join(os.path.dirname(here), "index.html")
open(out, "w", encoding="utf-8").write(html)
print(out, len(html.encode("utf-8")), "bytes")

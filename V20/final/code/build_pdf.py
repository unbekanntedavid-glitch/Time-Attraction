import re, markdown
from weasyprint import HTML
s = open('report_hu.md').read()
html = markdown.markdown(s, extensions=['tables', 'fenced_code', 'sane_lists'])
css = '''@page{size:A4;margin:14mm 12mm; @bottom-center{content: counter(page) " / " counter(pages); font-size:7pt}} body{font-family:"DejaVu Sans",sans-serif;font-size:8.6pt;line-height:1.33}
h1{font-size:15pt} h2{font-size:12.5pt;border-bottom:1px solid #999;margin-top:16px} h3{font-size:10.5pt}
table{border-collapse:collapse;width:100%;font-size:6.6pt;margin:6px 0} th,td{border:1px solid #bbb;padding:2px 3px;vertical-align:top;word-break:break-word}
th{background:#eee} code{font-size:7pt;word-break:break-all} li{margin:1px 0} img{width:100%;page-break-inside:avoid}'''
HTML(string=f'<html><head><meta charset="utf-8"><style>{css}</style></head><body>{html}</body></html>', base_url='.').write_pdf('report_hu.pdf')

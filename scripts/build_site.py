#!/usr/bin/env python3
"""Build a portable, offline static site from the canonical Markdown chapters."""
from __future__ import annotations
import html
from pathlib import Path
import re
from urllib.parse import urlsplit, urlunsplit

from bs4 import BeautifulSoup
import mistune

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
REPOSITORY_URL = 'https://github.com/FupaulHuang/drugseqpy-tutorials'
PAGES = [
    ('index', 'Learning path'), ('getting-started', 'Getting started'),
    ('input-preparation', 'Prepare inputs'), ('qc-workflow', 'Quality control'),
    ('mix-species-workflow', 'Mixed species'), ('multi-plate-comparison', 'Multi-plate comparison'),
    ('api-reference', 'API reference'), ('troubleshooting', 'Troubleshooting'),
]
ALIASES = {'qc_workflow_tutorial':'qc-workflow', 'mix_species_workflow_tutorial':'mix-species-workflow',
           'multi_plate_comparison_workflow_tutorial':'multi-plate-comparison', 'api_reference':'api-reference'}


def slug(text):
    return re.sub(r'-+', '-', re.sub(r'[^\w\s-]', '', text.lower()).replace(' ', '-')).strip('-')


def build():
    renderer = mistune.create_markdown(escape=False, plugins=['table', 'strikethrough', 'url'])
    for index, (name, label) in enumerate(PAGES):
        content = BeautifulSoup(renderer((DOCS / f'{name}.md').read_text()), 'html.parser')
        ids = {tag['id'] for tag in content.find_all(id=True)}
        for heading in content.find_all(re.compile('^h[1-4]$')):
            base = slug(heading.get_text(' ', strip=True))
            candidate, suffix = base, 1
            while candidate in ids:
                candidate = f'{base}-{suffix}'; suffix += 1
            heading['id'] = candidate; ids.add(candidate)
        for link in content.find_all('a', href=True):
            parts = urlsplit(link['href'])
            if not parts.scheme and not parts.netloc and parts.path.endswith('.md'):
                # Rewrite only local pages that are actually built.
                if Path(parts.path).stem in dict(PAGES) and '/' not in parts.path:
                    link['href'] = urlunsplit(parts._replace(path=parts.path[:-3]+'.html'))
        for img in content.find_all('img'):
            img['loading'] = 'lazy'; img['decoding'] = 'async'
        for table in content.find_all('table'):
            wrapper = content.new_tag('div', attrs={'class':'table-scroll','tabindex':'0','role':'region','aria-label':'Scrollable data table'})
            table.wrap(wrapper)
        toc = ''.join(f'<a href="#{h["id"]}">{html.escape(h.get_text(" ", strip=True))}</a>' for h in content.find_all('h2'))
        nav = ''.join(f'<a href="{key}.html"'+(' aria-current="page"' if key == name else '')+f'><span>{i+1:02}</span>{text}</a>' for i,(key,text) in enumerate(PAGES))
        neighbors = []
        if index:
            key,title = PAGES[index-1]; neighbors.append(f'<a href="{key}.html"><small>Previous</small>{title}</a>')
        if index < len(PAGES)-1:
            key,title=PAGES[index+1];neighbors.append(f'<a href="{key}.html"><small>Next</small>{title} →</a>')
        title = content.h1.get_text(' ',strip=True)
        page=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)} · drugseqpy tutorials</title>
<meta name="description" content="Connected drugseqpy tutorials for quality control, mixed-species analysis and multi-plate comparisons.">
<link rel="stylesheet" href="assets/site.css"><script src="assets/site.js" defer></script></head>
<body><a class="skip" href="#main">Skip to content</a>
<header class="topbar"><a class="brand" href="index.html"><span class="brand-mark">ds</span>drugseqpy <small>tutorials</small></a>
<div class="top-actions"><a href="{REPOSITORY_URL}/blob/main/docs/{name}.md">Markdown</a><a href="api-reference.html">API</a><button class="menu-toggle" aria-controls="chapter-nav" aria-expanded="false">Chapters</button></div></header>
<div class="layout"><aside class="sidebar" id="chapter-nav"><p class="nav-label">LEARN THE WORKFLOW</p><nav aria-label="Chapters">{nav}</nav>
<div class="sidebar-note">From counts to interpretation.<br>Sample 1–4 · Python workflows<br><a href="{REPOSITORY_URL}">GitHub repository ↗</a></div></aside>
<main id="main"><div class="chapter-label">DRUGSEQPY / {html.escape(label.upper())}</div><article>{content}</article>
<nav class="chapter-neighbors" aria-label="Previous and next chapters">{''.join(neighbors)}</nav>
<footer>drugseqpy tutorials · Worked examples and separate synthetic practice inputs.<br><a href="index.html">Learning path</a> · <a href="{name}.md">Read Markdown</a> · <a href="#main">Back to top ↑</a></footer></main>
<aside class="toc"><p class="nav-label">ON THIS PAGE</p><nav aria-label="On this page">{toc}</nav></aside></div></body></html>'''
        (DOCS / f'{name}.html').write_text(page)
    for old,new in ALIASES.items():
        (DOCS/f'{old}.html').write_text(f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta http-equiv="refresh" content="0;url={new}.html"><title>Chapter moved</title></head><body><p>Continue to <a href="{new}.html">{dict(PAGES)[new]}</a>.</p></body></html>''')
    (DOCS/'.nojekyll').touch()
    print(f'Built {len(PAGES)} chapters and {len(ALIASES)} compatibility redirects.')

if __name__ == '__main__':
    build()

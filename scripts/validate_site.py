#!/usr/bin/env python3
"""Validate local links, anchors, Python examples, privacy rules and synthetic files."""
import ast
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

from bs4 import BeautifulSoup
import mistune
from build_site import PAGES

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
errors = []
checked_links = 0
checked_python = 0
pages = {}
renderer = mistune.create_markdown(escape=False, plugins=['table', 'strikethrough', 'url'])


def fail(path, message):
    errors.append(f'{path.relative_to(ROOT)}: {message}')


def slug(text):
    return re.sub(r'-+', '-', re.sub(r'[^\w\s-]', '', text.lower()).replace(' ', '-')).strip('-')


def parse(path):
    if path in pages:
        return pages[path]
    text = path.read_text()
    soup = BeautifulSoup(renderer(text) if path.suffix == '.md' else text, 'html.parser')
    ids = [tag['id'] for tag in soup.find_all(id=True)]
    if path.suffix == '.md':
        for heading in soup.find_all(re.compile('^h[1-6]$')):
            candidate = slug(heading.get_text(' ', strip=True))
            original, n = candidate, 1
            while candidate in ids:
                candidate = f'{original}-{n}'; n += 1
            ids.append(candidate)
    else:
        duplicates = [key for key,n in Counter(ids).items() if n > 1]
        if duplicates:
            fail(path, f'duplicate anchors {duplicates}')
    pages[path] = soup, set(ids)
    return pages[path]

files = [ROOT/'README.md', ROOT/'VALIDATION.md'] + sorted((ROOT/'data').rglob('*.md')) + sorted(DOCS.glob('*.md')) + sorted(DOCS.glob('*.html'))
for path in files:
    if not path.exists():
        fail(path, 'missing document'); continue
    soup, ids = parse(path)
    for tag, attribute in [('a','href'), ('img','src'), ('script','src'), ('link','href')]:
        for item in soup.find_all(tag):
            target = item.get(attribute)
            if target is None:
                continue
            parts = urlsplit(target)
            if parts.scheme or parts.netloc:
                if tag in ('img','script','link'):
                    fail(path, f'nonlocal {tag} dependency')
                continue
            checked_links += 1
            resolved = (path.parent / unquote(parts.path)).resolve() if parts.path else path
            if not resolved.is_relative_to(ROOT):
                fail(path, f'link leaves repository: {target}');continue
            if not resolved.exists():
                fail(path, f'missing link target: {target}');continue
            if parts.fragment and resolved.suffix in ('.md','.html'):
                _, target_ids = parse(resolved)
                if unquote(parts.fragment) not in target_ids:
                    fail(path, f'missing anchor: {target}')
            if tag == 'img' and not item.get('alt'):
                fail(path, 'image has no alt text')
    if path.suffix == '.md' and path.name != 'api-reference.md':
        for number, code in enumerate(re.findall(r'^```python\n(.*?)^```', path.read_text(), re.M | re.S), 1):
            checked_python += 1
            try:
                ast.parse(code)
            except SyntaxError as error:
                fail(path, f'Python block {number}: {error}')
    if path.suffix == '.html' and path.stem in dict(PAGES):
        if len(soup.find_all('h1')) != 1:
            fail(path, 'expected one main heading')
        chapter_nav = soup.find('nav', attrs={'aria-label':'Chapters'})
        if not chapter_nav:
            fail(path, 'missing chapter navigation')
        elif [a.get('href') for a in chapter_nav.find_all('a')] != [f'{name}.html' for name, _ in PAGES]:
            fail(path, 'chapter navigation does not match the learning path')

# Scan publishable documents/assets metadata, never depend on private directories.
for path in list(DOCS.rglob('*')) + [ROOT/'README.md', ROOT/'VALIDATION.md']:
    if path.is_file() and path.suffix in ('.md','.html','.json','.css','.js','.svg'):
        text=path.read_text()
        if re.search(r'(?<![A-Za-z0-9])[Ff][1-4](?![A-Za-z0-9])', text):
            fail(path, 'unrevised sample identifier')
        if re.search(r'/(?:mnt|home|Users)/|nfs_share|file://|SOURCE_MANIFEST|SOURCE_PROVENANCE', text):
            fail(path, 'local source/provenance reference')
        if re.search(r'\b(?:GSE|GSM|SRR|SRP|ERR|ERP)\d+\b', text):
            fail(path, 'original dataset accession reference')
        if 'data:image/' in text:
            fail(path, 'embedded image instead of shared asset')

manifest=json.loads((ROOT/'data/practice-manifest.json').read_text())
if manifest['kind'] != 'synthetic':
    fail(ROOT/'data/practice-manifest.json','practice data not marked synthetic')
for name, expected in manifest['sha256'].items():
    path=ROOT/'data'/name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        fail(path,'practice checksum mismatch')
for path in (ROOT/'data').glob('sample_*/matrix/*.gz'):
    with gzip.open(path,'rb') as handle:
        handle.read()
bulk_manifest_path = ROOT/'data/bulk_practice/manifest.json'
bulk_manifest = json.loads(bulk_manifest_path.read_text())
if bulk_manifest.get('kind') != 'entirely synthetic':
    fail(bulk_manifest_path, 'bulk practice data not marked synthetic')
for name, expected in bulk_manifest['sha256'].items():
    path = bulk_manifest_path.parent/name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        fail(path, 'bulk practice checksum mismatch')
figures=list((DOCS/'assets/figures').glob('*.png'))
figure_manifest_path = DOCS/'assets/figures/manifest.json'
figure_manifest = json.loads(figure_manifest_path.read_text())
listed_assets = [entry['asset'] for entry in figure_manifest]
actual_assets = {p.relative_to(DOCS).as_posix() for p in figures}
if len(listed_assets) != len(set(listed_assets)) or set(listed_assets) != actual_assets:
    fail(figure_manifest_path, 'figure manifest and PNG assets differ or contain duplicates')
for chapter in {entry['chapter'] for entry in figure_manifest}:
    chapter_soup, _ = parse(DOCS/f'{chapter}.html')
    declared = [entry['asset'] for entry in figure_manifest if entry['chapter'] == chapter]
    rendered = [img.get('src') for img in chapter_soup.select('article img')]
    if declared != rendered:
        fail(figure_manifest_path, f'figure order differs from {chapter}.html')

if errors:
    print('\n'.join(errors))
    raise SystemExit(f'FAILED: {len(errors)} issues')
print(f'PASS: {len(files)} documents; {checked_links} local links/assets; {checked_python} Python examples; {len(figures)} figures; synthetic checksums.')

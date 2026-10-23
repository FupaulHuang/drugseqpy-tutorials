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

files = [ROOT/'README.md', ROOT/'VALIDATION.md', ROOT/'data/README.md'] + sorted(DOCS.glob('*.md')) + sorted(DOCS.glob('*.html'))
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
    if path.suffix == '.html' and path.name in ('index.html','qc-workflow.html','mix-species-workflow.html','multi-plate-comparison.html','api-reference.html'):
        if len(soup.find_all('h1')) != 1:
            fail(path, 'expected one main heading')
        if not soup.find('nav', attrs={'aria-label':'Chapters'}):
            fail(path, 'missing chapter navigation')

# Scan publishable documents/assets metadata, never depend on private directories.
for path in list(DOCS.rglob('*')) + [ROOT/'README.md', ROOT/'VALIDATION.md']:
    if path.is_file() and path.suffix in ('.md','.html','.json','.css','.js','.svg'):
        text=path.read_text()
        if re.search(r'(?<![A-Za-z0-9])[Ff][1-4](?![A-Za-z0-9])', text):
            fail(path, 'unrevised sample identifier')
        if re.search(r'/(?:mnt|home|Users)/|nfs_share|file://|SOURCE_MANIFEST|SOURCE_PROVENANCE', text):
            fail(path, 'local source/provenance reference')
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
figures=list((DOCS/'assets/figures').glob('*.png'))
if len(figures)!=68:
    fail(DOCS/'assets/figures', f'expected 68 figures, found {len(figures)}')

if errors:
    print('\n'.join(errors))
    raise SystemExit(f'FAILED: {len(errors)} issues')
print(f'PASS: {len(files)} documents; {checked_links} local links/assets; {checked_python} Python examples; {len(figures)} figures; synthetic checksums.')

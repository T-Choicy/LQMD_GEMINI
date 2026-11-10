"""Check the tracked release payload; regenerate hashes before tagging with --write-provenance."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = ROOT / 'PROVENANCE.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tracked_files():
    result = subprocess.run(['git', '-C', str(ROOT), 'ls-files', '-z'], capture_output=True, check=True)
    return sorted(p for p in result.stdout.decode('utf-8').split('\0') if p)


def main():
    paths = tracked_files()
    errors = []
    for name in paths:
        path = ROOT / name
        if re.search(r'(^|/)(tmp|temp|runs|data|x64|bin|__pycache__|\.venv|\.ai)(/|$)', name) or path.suffix.lower() in {'.exe', '.obj', '.mod', '.pdb', '.pyc', '.log', '.key', '.pem'} or path.name.startswith('.env'):
            errors.append(f'Forbidden payload: {name}')
        if path.stat().st_size > 10 * 1024**2:
            errors.append(f'Unexpected large payload (>10MiB project policy): {name}')
        if path.suffix.lower() == '.md':
            for match in re.finditer(r'\[[^\]]*\]\((?:<([^>]+)>|([^\)]+))\)', path.read_text(encoding='utf-8-sig')):
                link = match.group(1) or match.group(2)
                if not link or link.startswith('#') or urlsplit(link).scheme:
                    continue
                destination = unquote(link.split('#', 1)[0])
                if destination and not (path.parent / destination).exists():
                    errors.append(f'Broken local link: {name}: {link}')
        if name.startswith('fortran/') and path.suffix.lower() in {'.for', '.f'}:
            for line in path.read_text(encoding='utf-8-sig', errors='replace').splitlines():
                if line[:1].upper() in {'C', '*', '!'}:
                    continue
                include = re.match(r'\s*include\s+[\'\"]([^\'\"]+)', line, re.I)
                if include and not (path.parent / include.group(1)).exists():
                    errors.append(f'Missing include: {name}: {include.group(1)}')
    version = (ROOT / 'VERSION').read_text(encoding='utf-8-sig').strip()
    if not re.fullmatch(r'\d+\.\d+\.\d+', version):
        errors.append('Invalid VERSION')
    files = {name: digest(ROOT / name) for name in paths if name != 'PROVENANCE.json'}
    if '--write-provenance' in sys.argv:
        meta = json.loads(PROVENANCE.read_text(encoding='utf-8-sig')) if PROVENANCE.exists() else {}
        meta.update(version=version, files_sha256=files)
        PROVENANCE.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    else:
        meta = json.loads(PROVENANCE.read_text(encoding='utf-8-sig'))
        if meta.get('version') != version or meta.get('files_sha256') != files:
            errors.append('PROVENANCE differs from tracked payload; inspect changes then regenerate before release')
    if errors:
        print('\n'.join(errors), file=sys.stderr)
        return 1
    print(f'PASS v{version}: {len(paths)} tracked files, hashes, upload policy, local links and Fortran includes')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

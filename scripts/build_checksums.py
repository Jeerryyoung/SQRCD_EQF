from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "provenance" / "CHECKSUMS.sha256"
EXCLUDED_PARTS = {".git", ".venv", "__pycache__"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


files = []
for path in ROOT.rglob("*"):
    if not path.is_file() or path == OUT:
        continue
    if any(part in EXCLUDED_PARTS for part in path.parts):
        continue
    files.append(path)

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", encoding="utf-8", newline="\n") as stream:
    for path in sorted(files, key=lambda item: item.relative_to(ROOT).as_posix()):
        stream.write(f"{digest(path)}  {path.relative_to(ROOT).as_posix()}\n")

print(f"wrote {len(files)} checksums to {OUT}")

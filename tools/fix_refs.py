"""
fix_refs.py  —  PAPER_draft.md 참고문헌 정리 스크립트

변경 사항:
1. 본문 [19] → [1]  (ref 19 = Chen 2021이 OMP 설계 설명에 잘못 인용됨;
                       Thompson 2022 = ref 1이 올바른 인용)
2. 미인용 참고문헌 목록에서 삭제: 19, 29, 30, 31, 32, 33, 34, 35, 36, 37, 40, 42, 45
3. 나머지 참고문헌 연속 번호로 재정렬
"""

import re, sys, shutil, io
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

INFILE  = Path(r'D:\Claude-Code-R\LAMMPS-CO\PAPER_draft.md')
BACKUP  = INFILE.with_suffix('.md.bak')

# ── removal / mapping setup ──────────────────────────────────────────────────
REMOVE = {19, 29, 30, 31, 32, 33, 34, 35, 36, 37, 40, 42, 45}

# Build new numbering: keep refs not in REMOVE, renumber 1-N consecutively
kept = [i for i in range(1, 53) if i not in REMOVE]   # 39 refs
MAPPING = {old: new for new, old in enumerate(kept, start=1)}
# Special case: ref 19 in body text → ref 1 (Thompson 2022)
# (ref 19 is in REMOVE so it has no MAPPING entry; handle separately below)

# ── replacement function for body text ────────────────────────────────────
def replace_cit(m):
    nums = [int(n.strip()) for n in m.group(1).split(',')]
    # Sanity check: if none of the numbers are in 1-52, this isn't a citation
    if not any(1 <= n <= 52 for n in nums):
        return m.group(0)
    new_nums = []
    for n in nums:
        if n == 19:
            # Wrong citation in body → replace with ref 1 (Thompson 2022)
            new_nums.append(1)
        elif n in MAPPING:
            new_nums.append(MAPPING[n])
        elif n in REMOVE:
            pass  # drop this ref from the citation group (it's uncited anyway)
        else:
            new_nums.append(n)  # outside range — leave as-is
    if not new_nums:
        return m.group(0)   # nothing mapped — leave as-is (shouldn't happen)
    deduped = sorted(set(new_nums))
    return '[' + ','.join(str(x) for x in deduped) + ']'

CIT_RE = re.compile(r'(?<![A-Za-z0-9_])\[(\d+(?:,\s*\d+)*)\]')

# ── read file ────────────────────────────────────────────────────────────────
text = INFILE.read_text(encoding='utf-8')

# Split at ## References (find the LAST occurrence to be safe)
sep = '\n## References'
split_idx = text.rfind(sep)
if split_idx == -1:
    print("ERROR: '## References' section not found", file=sys.stderr)
    sys.exit(1)

body     = text[:split_idx]
refs_raw = text[split_idx:]

# ── process body: replace citations (skip code blocks) ──────────────────────
in_code  = False
new_body_lines = []
changes_body   = []

for lineno, line in enumerate(body.split('\n'), start=1):
    stripped = line.strip()
    if stripped.startswith('```'):
        in_code = not in_code
    if in_code or stripped.startswith('    ') and not stripped.startswith('    #'):
        # in code block — skip
        new_body_lines.append(line)
        continue
    new_line = CIT_RE.sub(replace_cit, line)
    if new_line != line:
        changes_body.append((lineno, line.strip(), new_line.strip()))
    new_body_lines.append(new_line)

new_body = '\n'.join(new_body_lines)

# ── process reference list ────────────────────────────────────────────────────
# Each ref entry: line matching  ^(\d+)\. (.+)$
# (single-line entries separated by blank lines)
REF_ENTRY = re.compile(r'^(\d+)\.\s+(.+)$')

new_refs_lines = []
refs_lines     = refs_raw.split('\n')
removed_refs   = []
kept_refs      = []
skip_blank     = False

i = 0
while i < len(refs_lines):
    line = refs_lines[i]
    m = REF_ENTRY.match(line)
    if m:
        old_num = int(m.group(1))
        content = m.group(2)
        if old_num in REMOVE:
            removed_refs.append(old_num)
            i += 1
            # skip immediately following blank line
            if i < len(refs_lines) and refs_lines[i].strip() == '':
                i += 1
            continue
        else:
            new_num = MAPPING.get(old_num, old_num)
            kept_refs.append((old_num, new_num))
            new_refs_lines.append(f'{new_num}. {content}')
    else:
        new_refs_lines.append(line)
    i += 1

new_refs = '\n'.join(new_refs_lines)

# ── report ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"Body text citation changes ({len(changes_body)} lines):")
for lineno, before, after in changes_body:
    print(f"  L{lineno}: {before}")
    print(f"       → {after}")

print(f"\nReference list: {len(removed_refs)} refs removed: {sorted(removed_refs)}")
print(f"Refs remaining: {len(kept_refs)} (total 39)")

# Show renumber summary for affected refs
print("\nRenumbering (old → new, changed only):")
for old, new in kept_refs:
    if old != new:
        print(f"  ref {old} → {new}")

# ── write ────────────────────────────────────────────────────────────────────
if '--dry-run' in sys.argv:
    print("\nDRY-RUN: file not modified.")
    sys.exit(0)

shutil.copy2(INFILE, BACKUP)
print(f"\nBackup saved: {BACKUP}")

INFILE.write_text(new_body + new_refs, encoding='utf-8')
print(f"Saved: {INFILE}")

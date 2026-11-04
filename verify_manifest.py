"""Check the frozen file inventory before running code that updates outputs."""
from pathlib import Path
import csv,hashlib
root=Path(__file__).resolve().parent
bad=[]
with (root/'MANIFEST.csv').open(encoding='utf-8-sig',newline='') as f:
 for row in csv.DictReader(f):
  p=root/row['path']
  if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:bad.append(row['path'])
if bad:raise SystemExit('MISSING OR CHANGED:\n'+'\n'.join(bad))
print('PASS: all files in MANIFEST.csv match; extra generated files are not checked')

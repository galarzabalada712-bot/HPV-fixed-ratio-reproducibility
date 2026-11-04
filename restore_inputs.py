"""Restore the two compressed frozen inputs and verify their original hashes."""
from pathlib import Path
import csv,gzip,hashlib,shutil
root=Path(__file__).resolve().parent
with (root/'reference/uncompressed_inputs_sha256.csv').open(encoding='utf-8',newline='') as f:
 for row in csv.DictReader(f):
  dest=root/row['path'];source=root/(row['path']+'.gz')
  if not dest.exists():
   with gzip.open(source,'rb') as a,dest.open('xb') as b:shutil.copyfileobj(a,b)
  digest=hashlib.sha256(dest.read_bytes()).hexdigest()
  if digest!=row['sha256']:raise SystemExit(f'Input hash mismatch: {row["path"]}; use a clean package copy.')
  print('PASS:',row['path'])

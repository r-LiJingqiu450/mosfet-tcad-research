"""Compatibility entry: explicit new S3 production freeze only."""
from pathlib import Path
import runpy
if __name__=='__main__':
    runpy.run_path(str(Path(__file__).resolve().parents[2]/'scripts/freeze_phase1.py'),run_name='__main__')

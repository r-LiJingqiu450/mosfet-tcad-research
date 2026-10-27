"""Release only completed, checked work; verify the ZIP itself after writing it."""
import argparse,hashlib,json,sys,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from verify_delivery import validate
from semiconductor_sim.common.jsonio import dumps


def source_digest(root):
    root=Path(root);digest=hashlib.sha256()
    paths=[]
    for dirname in ['semiconductor_sim','scripts']:
        paths.extend(p for p in (root/dirname).rglob('*.py') if '__pycache__' not in p.parts)
    paths.extend(root/n for n in ['pytest.ini','requirements.txt'])
    for p in sorted(paths):
        digest.update(p.relative_to(root).as_posix().encode());digest.update(p.read_bytes())
    return digest.hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('output');p.add_argument('--candidate',action='store_true');a=p.parse_args()
    validate(ROOT)
    clean=ROOT/'verification/clean_unpack_verification.json'
    if not a.candidate:
        if not clean.is_file():raise RuntimeError('Missing successful clean-unpack reproduction')
        record=json.loads(clean.read_text())
        if record.get('exit_code')!=0 or record.get('source_sha256')!=source_digest(ROOT):
            raise RuntimeError('Clean-unpack evidence does not match current source')
    paths=[]
    for f in ROOT.rglob('*'):
        if not f.is_file():continue
        rel=f.relative_to(ROOT);parts=rel.parts
        if any(x in parts for x in ['.cache','__pycache__','.pytest_cache']):continue
        include=(parts[0] in ['scripts','reports','baseline','history','verification']
             or rel.as_posix() in ['README.md','CHANGELOG.md','requirements.txt','environment_verified.txt',
                                 'pytest.ini','tests_phase1.txt','study_log.txt','changes.patch']
             or (parts[0]=='semiconductor_sim' and f.suffix=='.py')
             or rel.as_posix().startswith('semiconductor_sim/results/phase1/')
             or rel.as_posix()=='semiconductor_sim/results/golden/stage3.json')
        if include:paths.append((rel.as_posix(),f))
    entries=[dict(path=n,size_bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for n,f in sorted(paths)]
    manifest=dict(version='phase1-1.2',candidate=a.candidate,source_sha256=source_digest(ROOT),files=entries)
    target=Path(a.output).resolve();target.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
        for n,f in sorted(paths):z.write(f,'MOSFET_Phase1_v1_2/'+n)
        z.writestr('MOSFET_Phase1_v1_2/release_manifest.json',dumps(manifest,indent=2))
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(entries)+1
        for record in entries:
            data=z.read('MOSFET_Phase1_v1_2/'+record['path'])
            assert len(data)==record['size_bytes']
            assert hashlib.sha256(data).hexdigest()==record['sha256']
    print(dumps(dict(archive=str(target),files=len(entries)+1,size_bytes=target.stat().st_size,
                     sha256=hashlib.sha256(target.read_bytes()).hexdigest(),candidate=a.candidate)))

if __name__=='__main__':main()

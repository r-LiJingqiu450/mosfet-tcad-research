"""Run a command, sample its RSS every second, and count DEVSIM Newton 'Iteration' lines in its log at each sample."""
import subprocess, sys, time, os, json
log, out = sys.argv[1], sys.argv[2]
cmd = sys.argv[3:]
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
           DEVSIM_MATH_LIBS='/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3')
t0 = time.time()
with open(log, 'w') as f:
    p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, env=env, cwd='/home/claude/mosfet-tcad-research')
    samples = []
    while p.poll() is None:
        try:
            rss = int([l for l in open(f'/proc/{p.pid}/status') if l.startswith('VmRSS')][0].split()[1]) // 1024
        except Exception:
            rss = None
        txt = open(log, errors='ignore').read()
        samples.append(dict(t=round(time.time() - t0, 1), rss_MB=rss, iters=txt.count('Iteration:'), solves=txt.count('number of equations')))
        if time.time() - t0 > float(os.environ.get('MAXTIME_S', '1e9')):
            p.kill(); samples.append(dict(t=round(time.time() - t0, 1), note='stopped at MAXTIME_S')); break
        if rss is not None and rss > int(os.environ.get('MEMGUARD_MB', '5000')):
            p.kill(); samples.append(dict(t=round(time.time() - t0, 1), note='killed by memory guard')); break
        time.sleep(1.0)
json.dump(dict(cmd=cmd, exit=p.returncode, samples=samples), open(out, 'w'))
peak = max((s.get('rss_MB') or 0) for s in samples)
print('exit', p.returncode, 'wall', round(time.time() - t0), 's peak', peak, 'MB')
prev = None
for s in samples:
    if 'rss_MB' not in s:
        print(s); continue
    if prev is None or s['iters'] != prev['iters'] or (s['rss_MB'] or 0) - (prev['rss_MB'] or 0) > 150:
        print(s)
        prev = s

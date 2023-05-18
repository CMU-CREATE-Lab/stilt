#!/usr/bin/python3

import time

before = time.monotonic()
start_of_script = before

import glob, os, re, subprocess, sys, threading, traceback
import shutil

checkpoint = time.monotonic()
print(f"after loading non-pandas libs {(checkpoint - start_of_script)*1000:.0f} ms")

import pandas as pd

checkpoint = time.monotonic()
print(f"after loading pandas {(checkpoint - start_of_script)*1000:.0f} ms")

# Perf tests 3/3/23:
# 1 shard: 5.5 secs
# 2 shards: 3.0 secs
# 4 shards: 2.5 secs
# 8 shards: 1.8 secs
# 8 shards with parallel start: 1.3 secs.  HYSPLIT only took 0.8 secs
# 16 shards with parallel start: 1.16 secs.  HYSPLIT only took 0.575-0.697

# ThCall(func, *args, **kwargs) calls func(*args, **kwargs) in a separate thread
# value() waits for func to complete and returns its value
class ThCall(threading.Thread):
    def __init__(self, func, *args, **kwargs):
        self._exc_info = None
        self._output = {}
        def runner():
            try:
                retval = func(*args, **kwargs)
            except Exception as e:
                print(f'ThCall is relaying child exception {repr(e)} to parent', file=sys.stderr)
                sys.stderr.flush()
                self._output = {
                    "exception": e,
                    "traceback": traceback.format_exc()
                }
                return
            self._output = {"success": retval}
        super().__init__(target=runner)
        self.start()
    
    def value(self):
        if self.is_alive():
            self.join()
        if "exception" in self._output:
            e = self._output["exception"]
            print(f'ThCall is raising child exception in parent thread: {repr(e)}', file=sys.stderr)
            print(f'Child traceback: {self._output["traceback"]}', file=sys.stderr)
            sys.stderr.flush()
            raise e
        else:
            return self._output["success"]

n_shards = 16

# for file in ["CONC.CFG", "WARNING", "VMSDIST", "PARTICLE_STILT.DAT", "PARTICLE.DAT", "MESSAGE", "cdump"]:
#     try:
#         os.unlink(file)
#     except:
#         pass

# for dir in glob.glob("shard??"):
#     shutil.rmtree(dir, ignore_errors=True)

tmp_delete_dir = f"tmp_delete_{os.getpid()}_{threading.get_ident()}"
os.mkdir(tmp_delete_dir)
for path in ["CONC.CFG", "WARNING", "VMSDIST", "PARTICLE_STILT.DAT", "PARTICLE.DAT", "MESSAGE", "cdump"] + glob.glob("shard??"):
    if os.path.exists(path):
        os.rename(path, f"{tmp_delete_dir}/{path}")

ThCall(shutil.rmtree, tmp_delete_dir, ignore_errors=True)

#subprocess.check_output("rm -rf shard?? CONC.CFG WARNING VMSDIST PARTICLE_STILT.DAT PARTICLE.DAT MESSAGE cdump", shell=True)

checkpoint = time.monotonic()
print(f"after delete old files {(checkpoint - start_of_script)*1000:.0f} ms")



files_to_copy = [file for file in glob.glob("*") if not file.startswith("tmp_delete_")]
print("files_to_copy", files_to_copy)

checkpoint = time.monotonic()
print(f"after glob files {(checkpoint - start_of_script)*1000:.0f} ms")



# TODO: find MAXPAR and NUMPAR and divide them by # of shards.
setup_cfg = open("SETUP.CFG").read()
#print(f"setup_cfg is {setup_cfg}")
maxpar = int(re.search(r"^MAXPAR=(\d+)", setup_cfg, re.MULTILINE).group(1))
numpar = int(re.search(r"^NUMPAR=(\d+)", setup_cfg, re.MULTILINE).group(1))
#print(re.match(r"MAXPAR=(\d+)", "MAXPAR=3", re.MULTILINE))
shard_maxpar = int(maxpar/n_shards)
shard_numpar = int(numpar/n_shards)
print(f">>>      Changing MAXPAR from {maxpar} to {shard_maxpar}", file=sys.stderr)
print(f">>>      Changing NUMPAR from {numpar} to {shard_numpar}", file=sys.stderr)

shard_setup_cfg = re.sub(r"^MAXPAR=(\d+)", f"MAXPAR={shard_maxpar}", setup_cfg, flags=re.MULTILINE)
shard_setup_cfg = re.sub(r"^NUMPAR=(\d+)", f"NUMPAR={shard_numpar}", shard_setup_cfg, flags=re.MULTILINE)
assert(setup_cfg != shard_setup_cfg)

open("SETUP.CFG", "w").write(shard_setup_cfg)

checkpoint = time.monotonic()
print(f"after regex file editing and write new config {(checkpoint - start_of_script)*1000:.0f} ms")

# Set up shard dirs and kick off hysplits

def run_shard(shard_idx: int):
    before = time.monotonic()
    print(f"shard {shard_idx} starting at {(before - start_of_script)*1000:.0f} ms")
    shard_name = f"shard{shard_idx:02d}"
    os.mkdir(shard_name)
    for file_to_copy in files_to_copy:
        shutil.copy2(file_to_copy, shard_name)
    checkpoint = time.monotonic()
    print(f"shard {shard_idx} done copying at {(checkpoint - start_of_script)*1000:.0f} ms ({(checkpoint - before)*1000:.0f} ms elapsed)")
    # stdout = subprocess.check_output("./hycs_std", cwd=shard_name, stderr=subprocess.STDOUT)
    result = subprocess.run("./hycs_std", cwd=shard_name, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    checkpoint = time.monotonic()
    print(f"shard {shard_idx} done with subprocess_run at {(checkpoint - start_of_script)*1000:.0f} ms ({(checkpoint - before)*1000:.0f} ms elapsed)")
    if result.returncode != 0:
        msg = f"{result.returncode} with error: {result.stdout}"
        print(msg)
        raise RuntimeError(msg)
    print(f">>>      run_shard({shard_idx}) took {time.monotonic() - before:.3f} secs", file=sys.stderr)
    return result

# Start all shards
threadcalls = [ThCall(run_shard, shard_idx) for shard_idx in range(n_shards)]

particle_dfs = []

parsing_time = 0

for shard_idx in range(n_shards):
    stdout = threadcalls[shard_idx].value()
    parsing_time -= time.monotonic()

    print(f"%%% Shard {shard_idx} output {stdout}", file=sys.stderr)
    shard_name = f"shard{shard_idx:02d}"
    if shard_idx == 0:
        # Copy CONC.CFG from first shard, updating numpar and maxpar to total
        shard_conc_cfg = open(f"{shard_name}/CONC.CFG").read()
        #print(shard_conc_cfg)
        conc_cfg = re.sub(r"^ numpar\s*=\s*(\d+)", f" numpar = {shard_numpar * n_shards}", shard_conc_cfg, flags=re.MULTILINE)
        conc_cfg = re.sub(r"^ maxpar\s*=\s*(\d+)", f" maxpar = {shard_maxpar * n_shards}", conc_cfg, flags=re.MULTILINE)
        #print(conc_cfg)
        assert(conc_cfg != shard_conc_cfg)
        open("CONC.CFG", "w").write(conc_cfg)

    # TODO: synthesize from multiple particle files
    particle_stilt_header = open(f"{shard_name}/PARTICLE_STILT.DAT").readline()
    particle_df = pd.read_csv(f"{shard_name}/PARTICLE_STILT.DAT", sep=r'\s+', skiprows=1, header=None)
    #print(f"%%% Shard {shard_idx} has {len(particle_df)} rows", file=sys.stderr)
    particle_df[1] += shard_maxpar * shard_idx
    particle_dfs.append(particle_df)
    parsing_time += time.monotonic()

parsing_time -= time.monotonic()
particle_df = pd.concat(particle_dfs, axis=0)
# Sort descending according to time, ascending to particle number
particle_df.sort_values([0, 1], inplace=True, ascending=[False,True])

with open("PARTICLE_STILT.DAT", "w") as out:
    out.write(particle_stilt_header)
    particle_df.to_csv(out, sep=' ', index=None, header=None)

parsing_time += time.monotonic()
print(f">>>      parsing took {parsing_time:.3f} seconds", file=sys.stderr)

#print(f">>>      Collected has {len(particle_df)} rows", file=sys.stderr)
duration = time.monotonic() - before

print(f">>>      hycs_std_multi.py: computing with {n_shards} shards took {duration:.3f} seconds", file=sys.stderr)

#subprocess.check_output("cp PARTICLE_STILT.DAT CONC.CFG /tmp", shell=True)
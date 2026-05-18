#!/usr/bin/python3

import time

before = time.monotonic()
start_of_script = before

import glob, os, re, subprocess, sys, threading, traceback
import shutil

def log(message):
    # Prefix log messages with the current time, PID, and TID
    message = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {os.getpid()}.{threading.get_native_id()} server {message}"
    print(message)
    sys.stdout.flush()

def shard_message_tails(n_lines: int = 50) -> str:
    chunks = []
    for message_path in sorted(glob.glob("shard??/MESSAGE")):
        try:
            with open(message_path, "r", errors="replace") as f:
                tail = f.readlines()[-n_lines:]
        except Exception as e:
            chunks.append(f"===== {message_path} (failed to read: {repr(e)}) =====")
            continue

        chunks.append(f"===== {message_path} (last {n_lines} lines) =====")
        chunks.append("".join(tail).rstrip("\n"))

    return "\n".join(chunks)

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

n_shards = 26

tmp_delete_dir = f"tmp_delete_{os.getpid()}_{threading.get_ident()}"
os.mkdir(tmp_delete_dir)
for path in ["CONC.CFG", "WARNING", "VMSDIST", "PARTICLE_STILT.DAT", "PARTICLE.DAT", "MESSAGE", "cdump"] + glob.glob("shard??"):
    if os.path.exists(path):
        os.rename(path, f"{tmp_delete_dir}/{path}")

ThCall(shutil.rmtree, tmp_delete_dir, ignore_errors=True)

files_to_link = [file for file in glob.glob("*") if not file.startswith("tmp_delete_")]

setup_cfg = open("SETUP.CFG").read()
maxpar = int(re.search(r"^MAXPAR=(\d+)", setup_cfg, re.MULTILINE).group(1))
numpar = int(re.search(r"^NUMPAR=(\d+)", setup_cfg, re.MULTILINE).group(1))
shard_maxpar = int(maxpar/n_shards)
shard_numpar = int(numpar/n_shards)

shard_setup_cfg = re.sub(r"^MAXPAR=(\d+)", f"MAXPAR={shard_maxpar}", setup_cfg, flags=re.MULTILINE)
shard_setup_cfg = re.sub(r"^NUMPAR=(\d+)", f"NUMPAR={shard_numpar}", shard_setup_cfg, flags=re.MULTILINE)
assert(setup_cfg != shard_setup_cfg)

open("SETUP.CFG", "w").write(shard_setup_cfg)

# Set up shard dirs and kick off hysplits
def run_shard(shard_idx: int):
    before = time.monotonic()
    shard_name = f"shard{shard_idx:02d}"
    os.mkdir(shard_name)
    for file_to_link in files_to_link:
        source_path = os.path.join(os.getcwd(), file_to_link)
        target_path = os.path.join(shard_name, file_to_link)
        try:
            os.symlink(source_path, target_path)
        except OSError:
            shutil.copy2(source_path, target_path)
    exe_path = os.path.join(os.getcwd(), shard_name, "hycs_std")
    try:
        result = subprocess.run("./hycs_std", cwd=shard_name, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except Exception as e:
        log(f"{exe_path} failed with exception {repr(e)}")
        raise
    if result.returncode != 0:
        msg = f"{exe_path} failed with status {result.returncode} and error {result.stdout}"
        log(msg)
        raise RuntimeError(msg)
    log(f"run_shard({shard_idx}) took {time.monotonic() - before:.3f} secs")
    return result

# Start all shards
threadcalls = [ThCall(run_shard, shard_idx) for shard_idx in range(n_shards)]

merge_cmdline = [
    f"{os.path.dirname(__file__)}/merge_particle_stilt_files",
    "PARTICLE_STILT.DAT"
]

for shard_idx in range(n_shards):
    stdout = threadcalls[shard_idx].value()

    shard_name = f"shard{shard_idx:02d}"
    if shard_idx == 0:
        # Copy CONC.CFG from first shard, updating numpar and maxpar to total
        shard_conc_cfg = open(f"{shard_name}/CONC.CFG").read()
        conc_cfg = re.sub(r"^ numpar\s*=\s*(\d+)", f" numpar = {shard_numpar * n_shards}", shard_conc_cfg, flags=re.MULTILINE)
        conc_cfg = re.sub(r"^ maxpar\s*=\s*(\d+)", f" maxpar = {shard_maxpar * n_shards}", conc_cfg, flags=re.MULTILINE)
        assert(conc_cfg != shard_conc_cfg)
        open("CONC.CFG", "w").write(conc_cfg)

    merge_cmdline.append(f"{shard_name}/PARTICLE_STILT.DAT")
    merge_cmdline.append(f"{shard_maxpar * shard_idx}")

try:
    subprocess.check_output(merge_cmdline, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
    output = e.output.decode("utf-8", errors="replace") if e.output else ""
    shard_tails = shard_message_tails(50)
    msg = f"merge_particle_stilt_files failed with status {e.returncode}. Output:\n{output}"
    if shard_tails:
        msg = f"{msg}\n\nShard MESSAGE tails:\n{shard_tails}"
    log(msg)
    raise RuntimeError(msg)

if not os.path.exists("PARTICLE_STILT.DAT") or os.path.getsize("PARTICLE_STILT.DAT") == 0:
    msg = "merge_particle_stilt_files completed but PARTICLE_STILT.DAT was not created or is empty"
    shard_tails = shard_message_tails(50)
    if shard_tails:
        msg = f"{msg}\n\nShard MESSAGE tails:\n{shard_tails}"
    log(msg)
    raise RuntimeError(msg)

duration = time.monotonic() - before

print(f">>>      hycs_std_multi.py: computing with {n_shards} shards took {duration:.3f} seconds", file=sys.stderr)

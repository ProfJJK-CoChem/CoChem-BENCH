import sys
import os
from pathlib import Path

repo_path = r"D:\Gdrive\__CoChem\GitHub-Repo"
if repo_path not in sys.path:
    sys.path.insert(0, repo_path)

import time
import numpy as np

# Adjust the import to match the actual path structure. 
# The folder is "CoChem-KINETIC", which Python cannot import directly with a hyphen.
# Wait, "CoChem-KINETIC\kinetic_core" means the package might be "kinetic_core" or we append "CoChem-KINETIC" to sys.path.
sys.path.insert(0, os.path.join(repo_path, "CoChem-KINETIC"))
from kinetic_core.cochem_pes_store import PESStore

def test_write_speed():
    store_path = Path("test_1gb.h5")
    if store_path.exists():
        store_path.unlink()
    
    store = PESStore(store_path)
    
    n_steps = 10000
    N_ATOMS = 1500  # 1500 * 3 * 8 * 10000 = ~360 MB for coords, ~360 MB for grad. Total ~720MB. 
    # Let's make it 2000 atoms to hit close to 1GB.
    N_ATOMS = 2200
    
    print("Generating 1GB dummy trajectory dataset in RAM...")
    coords = np.random.rand(n_steps, N_ATOMS, 3)
    energy = np.random.rand(n_steps)
    gradient = np.random.rand(n_steps, N_ATOMS, 3)
    
    size_mb = (coords.nbytes + energy.nbytes + gradient.nbytes) / (1024 * 1024)
    print(f"Data size: {size_mb:.2f} MB")
    
    print("Writing to PESStore...")
    start_time = time.time()
    store.append_batch(coords, energy, gradient_batch=gradient)
    end_time = time.time()
    
    write_time = end_time - start_time
    speed_mb_s = size_mb / write_time
    
    print(f"Write time: {write_time:.2f} s")
    print(f"Write speed: {speed_mb_s:.2f} MB/s")
    
    if store_path.exists():
        store_path.unlink()
        
    if speed_mb_s < 100:
        print("FAIL: Write speed < 100 MB/s")
        sys.exit(1)
    else:
        print("PASS: Write speed > 100 MB/s")
        sys.exit(0)

if __name__ == "__main__":
    test_write_speed()

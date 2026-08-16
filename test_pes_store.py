import os
import numpy as np
from cochem_bench.core.storage_io import PESStore

def test():
    print("Testing PESStore...")
    filepath = "test_pes.h5"
    if os.path.exists(filepath):
        os.remove(filepath)
        
    store = PESStore(filepath)
    coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    energies = np.array([-1.0, -1.2])
    attrs = {"method": "DFT", "basis": "6-31G*"}
    
    print("Saving PES...")
    store.save_pes("mol1", coords, energies, attrs)
    
    print("Loading PES...")
    loaded_coords, loaded_energies, loaded_attrs = store.load_pes("mol1")
    
    assert np.allclose(coords, loaded_coords), "Coordinates mismatch"
    assert np.allclose(energies, loaded_energies), "Energies mismatch"
    assert attrs == loaded_attrs, "Attributes mismatch"
    
    print("Test passed successfully. Import and execution SUCCESS.")
    
    if os.path.exists(filepath):
        os.remove(filepath)

if __name__ == "__main__":
    test()

import os
import sys

# Ensure cochem_bench is in path
sys.path.insert(0, os.path.abspath('D:\\Gdrive\\__CoChem\\GitHub-Repo\\CoChem-BENCH'))

try:
    from cochem_bench.hardware import setup_parallelism, HPCConfig
    
    print("Importing cochem_bench.hardware successful!")

    # Test default
    config = setup_parallelism()
    print(f"Default config: world_size={config.world_size}, rank={config.rank}, local_rank={config.local_rank}, is_master={config.is_master}")

    # Set MPI environment variables and test again
    os.environ["OMPI_COMM_WORLD_SIZE"] = "4"
    os.environ["OMPI_COMM_WORLD_RANK"] = "2"
    os.environ["OMPI_COMM_WORLD_LOCAL_RANK"] = "1"
    
    # Re-instantiate directly to bypass singleton
    config_mpi = HPCConfig()
    print(f"MPI config: world_size={config_mpi.world_size}, rank={config_mpi.rank}, local_rank={config_mpi.local_rank}, is_master={config_mpi.is_master}")
    print(f"Node count (2 ranks per node): {config_mpi.get_node_count(2)}")

    print("Functional execution test passed!")

except Exception as e:
    print(f"Failed: {e}")
    sys.exit(1)

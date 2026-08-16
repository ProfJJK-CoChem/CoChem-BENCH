import time
import sys

# Attempt to import the target module
try:
    from cochem_bench.intake.caching import get_pubchem_records
except ImportError as e:
    print(f"ImportError: {e}")
    sys.exit(1)

def test_cache_miss_bomb():
    # 100 CIDs
    cids = list(range(2244, 2344))
    
    # Pass 1: Cache miss
    start_time_1 = time.perf_counter()
    pass1_results = get_pubchem_records(cids)
    end_time_1 = time.perf_counter()
    pass1_duration = end_time_1 - start_time_1
    
    print(f"Pass 1 (Network) Duration: {pass1_duration:.4f} seconds")
    
    # Pass 2: Cache hit
    start_time_2 = time.perf_counter()
    pass2_results = get_pubchem_records(cids)
    end_time_2 = time.perf_counter()
    pass2_duration = end_time_2 - start_time_2
    
    print(f"Pass 2 (Cache) Duration: {pass2_duration:.4f} seconds")
    
    # Assertions
    assert pass2_duration < 1.0, f"Cache retrieval took > 1.0 seconds: {pass2_duration:.4f} s"
    print("Latency Assertion Passed: Cache retrieval is < 1.0 seconds.")
    
    speedup = pass1_duration / pass2_duration
    print(f"Speedup: {speedup:.2f}x")
    assert pass2_duration * 100 <= pass1_duration, f"Pass 2 is not at least 100x faster than pass 1! Speedup: {speedup:.2f}x"
    print("Validation Passed: Pass 2 is at least 100x faster than Pass 1.")

if __name__ == "__main__":
    test_cache_miss_bomb()

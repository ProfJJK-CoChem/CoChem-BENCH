import time
import sys

# Attempt to import the target module
try:
    from cochem_base.core.swarm_messaging import generate_teamwork_preview
except ImportError as e:
    print(f"ImportError: {e}")
    sys.exit(1)

def test_swarm_latency_bomb():
    print("Running 100-Agent Bomb test...")
    start_time = time.perf_counter()
    
    # Run the teamwork preview
    result = generate_teamwork_preview(100)
    
    time_elapsed = time.perf_counter() - start_time
    print(f"Time elapsed: {time_elapsed:.4f} seconds")
    
    assert time_elapsed < 5.0, "Latency exceeded 5 seconds!"
    assert "Task_List.md" in result, "Output does not contain Task_List.md"
    assert "Task 100" in result, "Output does not contain 100 tasks"
    print("Validation Passed: Latency < 5.0s and output generated correctly.")

if __name__ == "__main__":
    test_swarm_latency_bomb()

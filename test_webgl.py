import sys

try:
    from cochem_bench.interfaces import WebGLPerformanceTracker, benchmark_fps
    print("SUCCESS: Successfully imported WebGLPerformanceTracker and benchmark_fps.")
    
    tracker = WebGLPerformanceTracker(fps_threshold=30.0)
    print(f"Initial mode: {tracker.get_rendering_mode()}")
    
    # Run a short evaluation that will likely pass WebGL
    print("Evaluating performance (low polygons)...")
    is_svg, fps = tracker.evaluate_performance(benchmark_duration=0.2, num_polygons=1)
    print(f"FPS: {fps:.2f}, Fallback to SVG: {is_svg}, Mode: {tracker.get_rendering_mode()}")
    
    # Run a short evaluation that might trigger SVG
    print("Evaluating performance (high polygons)...")
    is_svg, fps = tracker.evaluate_performance(benchmark_duration=0.2, num_polygons=1000000)
    print(f"FPS: {fps:.2f}, Fallback to SVG: {is_svg}, Mode: {tracker.get_rendering_mode()}")
    
    print("FUNCTIONAL TEST PASSED")
    sys.exit(0)
    
except Exception as e:
    print(f"FAILED: {e}")
    sys.exit(1)

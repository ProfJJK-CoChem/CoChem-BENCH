import logging
logging.basicConfig(level=logging.INFO)
import sys

try:
    from cochem_bench.interfaces import WebGLPerformanceTracker, benchmark_fps
    logging.info("SUCCESS: Successfully imported WebGLPerformanceTracker and benchmark_fps.")
    
    tracker = WebGLPerformanceTracker(fps_threshold=30.0)
    logging.info(f"Initial mode: {tracker.get_rendering_mode()}")
    
    # Run a short evaluation that will likely pass WebGL
    logging.info("Evaluating performance (low polygons)...")
    is_svg, fps = tracker.evaluate_performance(benchmark_duration=0.2, num_polygons=1)
    logging.info(f"FPS: {fps:.2f}, Fallback to SVG: {is_svg}, Mode: {tracker.get_rendering_mode()}")
    
    # Run a short evaluation that might trigger SVG
    logging.info("Evaluating performance (high polygons)...")
    is_svg, fps = tracker.evaluate_performance(benchmark_duration=0.2, num_polygons=1000000)
    logging.info(f"FPS: {fps:.2f}, Fallback to SVG: {is_svg}, Mode: {tracker.get_rendering_mode()}")
    
    logging.info("FUNCTIONAL TEST PASSED")
    sys.exit(0)
    
except Exception as e:
    logging.info(f"FAILED: {e}")
    sys.exit(1)

import logging
logging.basicConfig(level=logging.INFO)
import sys
import time

def test_webgl_fallback():
    logging.info("Running WebGL Headless Rendering FPS Audit...")
    try:
        from cochem_bench.interfaces import webgl_perf
    except ImportError as e:
        logging.info(f"[FAIL] ImportError: {e}")
        sys.exit(1)

    logging.info("Feeding 50,000-atom protein .pdb to the Web-MCP interface...")
    try:
        generator = webgl_perf.UI_Generator()
        
        logging.info("Forcing renderer into SVG mode...")
        generator.set_mode('SVG')
        
        # Simulate rendering and poll latency
        generator.render("50k_atom_protein.pdb")
        fps = generator.get_fps()
        logging.info(f"SVG Mode FPS: {fps}")
        
        if fps < 1.0:
            logging.info("Rendering latency detected. Catastrophically stalled (< 1 FPS).")
        
        # Check if fallback sequence was triggered and WebGL mode engaged
        current_mode = generator.get_mode()
        logging.info(f"Current UI Mode: {current_mode}")
        
        assert current_mode == 'WebGL', "Fallback to WebGL was not triggered."
        logging.info("[PASS] Fallback sequence executed and WebGL mode is engaged.")
        
        logging.info("Rendering WebGL topology to benchmark array buffer performance...")
        generator.render("50k_atom_protein.pdb")
        webgl_fps = generator.get_fps()
        logging.info(f"WebGL Mode FPS: {webgl_fps}")
        
        assert webgl_fps > 30.0, f"WebGL failed to maintain > 30 FPS! Actual: {webgl_fps}"
        logging.info("[PASS] WebGL hardware acceleration maintained >30 FPS constraint.")
        
    except Exception as e:
        logging.info(f"[FAIL] Execution error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    test_webgl_fallback()

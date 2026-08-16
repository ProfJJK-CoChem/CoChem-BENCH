import sys
import time

def test_webgl_fallback():
    print("Running WebGL Headless Rendering FPS Audit...")
    try:
        from cochem_bench.interfaces import webgl_perf
    except ImportError as e:
        print(f"[FAIL] ImportError: {e}")
        sys.exit(1)

    print("Feeding 50,000-atom protein .pdb to the Web-MCP interface...")
    try:
        generator = webgl_perf.UI_Generator()
        
        print("Forcing renderer into SVG mode...")
        generator.set_mode('SVG')
        
        # Simulate rendering and poll latency
        generator.render("50k_atom_protein.pdb")
        fps = generator.get_fps()
        print(f"SVG Mode FPS: {fps}")
        
        if fps < 1.0:
            print("Rendering latency detected. Catastrophically stalled (< 1 FPS).")
        
        # Check if fallback sequence was triggered and WebGL mode engaged
        current_mode = generator.get_mode()
        print(f"Current UI Mode: {current_mode}")
        
        assert current_mode == 'WebGL', "Fallback to WebGL was not triggered."
        print("[PASS] Fallback sequence executed and WebGL mode is engaged.")
        
        print("Rendering WebGL topology to benchmark array buffer performance...")
        generator.render("50k_atom_protein.pdb")
        webgl_fps = generator.get_fps()
        print(f"WebGL Mode FPS: {webgl_fps}")
        
        assert webgl_fps > 30.0, f"WebGL failed to maintain > 30 FPS! Actual: {webgl_fps}"
        print("[PASS] WebGL hardware acceleration maintained >30 FPS constraint.")
        
    except Exception as e:
        print(f"[FAIL] Execution error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    test_webgl_fallback()

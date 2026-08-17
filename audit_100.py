import logging
logging.basicConfig(level=logging.INFO)
import sys
import traceback

def run_audit():
    try:
        # The prompt specifies Target: cochem_bench.literature.openalex_graph
        try:
            from cochem_bench.literature import openalex_graph
        except ImportError:
            try:
                from cochem_bench.literature import openalex_graph
            except ImportError as e:
                logging.info("[FAIL] The target module 'cochem_bench.literature.openalex_graph' (or 'cochem_bench.literature.openalex_graph') could not be imported.")
                logging.info(f"Details: {e}")
                sys.exit(1)

        # Ensure GraphDepthWarning exists
        if not hasattr(openalex_graph, "GraphDepthWarning"):
            logging.info("[FAIL] GraphDepthWarning is not defined in openalex_graph.")
            sys.exit(1)

        import warnings
        
        logging.info("Invoking openalex_graph with depth=4...")
        try:
            # We are guessing the builder's name, e.g. build_graph or fetch_network
            # We'll just inspect the module
            func = getattr(openalex_graph, "build_graph", None)
            if not func:
                logging.info("[FAIL] Could not find the graph builder function in openalex_graph.")
                sys.exit(1)
            
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                graph = func("W2018541907", depth=4)
                
                # Check for GraphDepthWarning
                warning_found = any(issubclass(warn.category, openalex_graph.GraphDepthWarning) for warn in w)
                
                if not warning_found:
                    logging.info("[FAIL] depth=4 did not issue GraphDepthWarning.")
                    sys.exit(1)
                    
                if not graph or not graph.get('nodes'):
                    logging.info("[FAIL] Graph is empty or not returned.")
                    sys.exit(1)
                    
                if not graph.get('edges'):
                    logging.info("[FAIL] Graph has no edges, meaning no references were fetched.")
                    sys.exit(1)
                    
                if "GraphDepthWarning" not in graph.get('metadata', {}).get('warnings', []):
                    logging.info("[FAIL] GraphDepthWarning not appended to metadata.")
                    sys.exit(1)
                    
            logging.info("[PASS] GraphDepthWarning was issued and graph was truncated correctly.")
            sys.exit(0)
        except Exception as e:
            logging.info(f"[FAIL] Unexpected error: {e}")
            traceback.print_exc()
            sys.exit(1)

    except Exception as e:
        logging.info(f"[FAIL] Audit script failed to run properly: {e}")
        sys.exit(1)

if __name__ == "__main__":
    run_audit()

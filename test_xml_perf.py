import logging
logging.basicConfig(level=logging.INFO)
import sys
import os

# Add CoChem-BENCH to sys.path so we can import cochem_bench
sys.path.insert(0, r"D:\Gdrive\__CoChem\GitHub-Repo\CoChem-BENCH")

from cochem_bench.literature.xml_performance import XMLParsingBenchmark, validate_xml_memory_usage

def main():
    logging.info("Testing import of XMLParsingBenchmark...")
    
    # Create a dummy XML file
    dummy_xml_path = "dummy_test.xml"
    with open(dummy_xml_path, "w") as f:
        f.write("<root>")
        for i in range(100):
            f.write(f"<item id='{i}'>Item {i}</item>")
        f.write("</root>")
        
    try:
        logging.info("Testing validate_xml_memory_usage...")
        is_valid, stats = validate_xml_memory_usage(dummy_xml_path, memory_limit_mb=10.0)
        logging.info(f"Validation successful: {is_valid}")
        logging.info(f"Stats: {stats}")
        
        logging.info("Testing memory limit exceeded simulation...")
        # Very small limit to force failure
        is_valid_fail, stats_fail = validate_xml_memory_usage(dummy_xml_path, memory_limit_mb=0.000001)
        logging.info(f"Validation (expected fail): {is_valid_fail}")
        logging.info(f"Stats (fail): {stats_fail}")
        
    finally:
        if os.path.exists(dummy_xml_path):
            os.remove(dummy_xml_path)

if __name__ == "__main__":
    main()

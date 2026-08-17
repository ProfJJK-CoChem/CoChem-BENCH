import logging
logging.basicConfig(level=logging.INFO)
import os
import sys

logging.info("--- Running CoChem-BENCH Adversarial Audit 93 ---")

# Add CoChem-NODE and CoChem-BENCH to path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_root, 'CoChem-NODE'))
sys.path.append(os.path.join(project_root, 'CoChem-BENCH'))

def run_audit():
    try:
        from cochem_slurm_templater import SlurmTemplater
        
        # Mocking a configuration where node max cores is 64
        class MockHardware:
            cpu_cores = 64
            ram_mb = 128000
            
        class MockConfig:
            hardware = MockHardware()
            walltime_budgets = {"T1-30min": "00:30:00"}
            
        templater = SlurmTemplater(config=MockConfig())
        
        script = templater.render_job(
            job_name="Audit93_Job",
            work_dir=".",
            execution_command="orca input.inp",
            requested_cores=256
        )
        
        logging.info("\n--- Generated Batch Script ---")
        logging.info(script)
        logging.info("------------------------------")
        
        if "mpirun -np 256" in script or "srun" in script:
            logging.info("\n[VERIFICATION] PASS: Multi-Node HPC Routing identified, mpirun wrapper present.")
        elif "OMP_NUM_THREADS=256" in script:
            logging.info("\n[VERIFICATION] FAIL: System attempted to spawn 256 OpenMP threads. Severe context-switching overhead!")
        else:
            logging.info("\n[VERIFICATION] FAIL: Neither mpirun -np 256 nor OMP_NUM_THREADS=256 found. Check resource throttling.")
            if "OMP_NUM_THREADS=64" in script:
                logging.info("Note: The system throttled the request down to 64 cores instead of deploying a Multi-Node job.")

    except Exception as e:
        logging.info(f"\n[ERROR] Audit execution failed: {e}")

if __name__ == "__main__":
    run_audit()

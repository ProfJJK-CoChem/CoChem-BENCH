import streamlit as st
import subprocess
import os
import sys
import psutil
import atexit
import hashlib
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)

st.set_page_config(page_title="CoChem-BENCH - Native Pipeline UI", layout="wide")

def kill_zombie_processes() -> None:
    target_procs = ['orca', 'xtb', 'mpi', 'crest']
    for proc in psutil.process_iter(['name']):
        try:
            name = proc.info.get('name')
            if name and any(target in name.lower() for target in target_procs):
                proc.terminate()
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess) as e:
            logging.warning(f"Failed to terminate zombie process: {e}")
            continue
atexit.register(kill_zombie_processes)

# Initialize Session State
if "run_history" not in st.session_state:
    st.session_state.run_history = []
if "latest_stdout" not in st.session_state:
    st.session_state.latest_stdout = ""

st.title("🔬 CoChem-BENCH Control Panel")
st.markdown("This UI executes raw, heavy mathematical payloads natively.")

tab1, tab2, tab3 = st.tabs(["⚙️ Configuration", "📜 Logs", "📊 Results"])

with tab1:
    with st.sidebar:
        st.header("Pipeline Configuration")
        target_smiles = st.text_input("Target SMILES", "CCO", help="Input a valid SMILES string for the target molecule.")
        run_mode = st.selectbox("Execution Mode", ["Fast", "Accurate"], help="Select the fidelity of the quantum chemistry methods.")
        
        st.markdown("---")
        st.markdown("### Process Management")
        if st.button("🧹 Force Cleanup Zombies"):
            kill_zombie_processes()
            st.success("Zombie processes terminated.")

    if st.button("🚀 Execute Default Pipeline", use_container_width=True):
        with st.spinner(f"Triggering quantum physics executor for {target_smiles}..."):
            st.info("Initiating Physical Math Execution Pipeline...")
            
            module_dir = Path(__file__).resolve().parent
            dispatcher_path = module_dir / 'cochem_bench' / 'dispatcher.py'
            
            if not dispatcher_path.exists():
                st.error(f"Dispatcher not found at {dispatcher_path}")
            else:
                env = os.environ.copy()
                env["COCHEM_TARGET_H5"] = os.path.join(os.getcwd(), "landscape.h5")
                env["COCHEM_TARGET_SMILES"] = target_smiles
                env["COCHEM_RUN_MODE"] = run_mode
                
                try:
                    cmd = [sys.executable, str(dispatcher_path)]
                    result = subprocess.run(
                        cmd, 
                        capture_output=True, 
                        text=True, 
                        check=True, 
                        timeout=3600, 
                        cwd=str(module_dir),
                        env=env
                    )
                    
                    st.session_state.latest_stdout = result.stdout
                    st.session_state.run_history.append({"smiles": target_smiles, "mode": run_mode, "status": "Success"})
                    
                    st.success(f"✅ Execution Completed Natively for {target_smiles}. CPU load generated.")
                    
                    # [ANTI-SPOOFING] Removed fake generation of physical_output.out.
                    # The dispatcher must produce the actual files.
                        
                except subprocess.TimeoutExpired:
                    st.error("Execution timed out. Purging zombies.")
                    st.session_state.run_history.append({"smiles": target_smiles, "mode": run_mode, "status": "Timeout"})
                    kill_zombie_processes()
                except subprocess.CalledProcessError as e:
                    st.warning(f"Execution finished with non-zero exit code: {e.returncode}")
                    st.session_state.latest_stdout = e.stdout + "\n" + e.stderr
                    st.session_state.run_history.append({"smiles": target_smiles, "mode": run_mode, "status": f"Failed ({e.returncode})"})
                    kill_zombie_processes()
                except Exception as e:
                    st.error(f"Pipeline crashed during physical execution: {str(e)}")
                    st.session_state.run_history.append({"smiles": target_smiles, "mode": run_mode, "status": "Error"})
                    kill_zombie_processes()

with tab2:
    st.subheader("Execution Logs")
    if st.session_state.latest_stdout:
        st.code(st.session_state.latest_stdout, language="text")
    else:
        st.info("No logs available yet. Run a pipeline first.")

with tab3:
    st.subheader("Run History")
    if st.session_state.run_history:
        st.table(st.session_state.run_history)
    else:
        st.info("No runs have been executed in this session.")

import logging
logging.basicConfig(level=logging.INFO)
import sys
import base64
import os
import tiktoken

def count_tokens(text: str, model: str = "gpt-4") -> int:
    try:
        encoding = tiktoken.encoding_for_model(model)
    except KeyError:
        encoding = tiktoken.get_encoding("cl100k_base")
    return len(encoding.encode(text))

def run_audit():
    try:
        from cochem_bench import llm_context
    except ImportError:
        logging.info("[FAIL] cochem_bench.llm_context module not found.")
        sys.exit(1)

    logging.info("Generating 40,000+ tokens of random Base64 strings...")
    system_prompt = "[SYSTEM PROMPT] You are a helpful assistant."
    last_user_request = "[LAST USER REQUEST] Please fix the CFOUR syntax."
    
    middle_logs = []
    # Make sure we generate well over 40k tokens
    for _ in range(2000):
        # 10 words per log
        log_entry = " ".join([base64.b64encode(os.urandom(6)).decode('utf-8') for _ in range(10)])
        middle_logs.append(log_entry)
        
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": "Initial request"},
    ]
    for log in middle_logs:
        messages.append({"role": "assistant", "content": "Log dump"})
        messages.append({"role": "user", "content": log})
        
    messages.append({"role": "user", "content": last_user_request})
    
    total_input_tokens = sum(count_tokens(m["content"]) for m in messages)
    logging.info(f"Total tokens before processing: {total_input_tokens}")
    
    try:
        result_messages = llm_context.build_message_queue(messages)
    except Exception as e:
        logging.info(f"[FAIL] Error running build_message_queue: {e}")
        sys.exit(1)
        
    total_tokens = sum(count_tokens(m["content"]) for m in result_messages)
    logging.info(f"Resulting tokens after processing: {total_tokens}")
    
    if total_tokens > 32000:
        logging.info("[FAIL] Resulting token array > 32000.")
        sys.exit(1)
        
    if result_messages[0]["content"] != system_prompt:
        logging.info("[FAIL] System prompt not preserved.")
        sys.exit(1)
        
    if result_messages[-1]["content"] != last_user_request:
        logging.info("[FAIL] Last user request not preserved.")
        sys.exit(1)
        
    logging.info("[PASS] LLM context degradation audit passed.")

if __name__ == "__main__":
    run_audit()

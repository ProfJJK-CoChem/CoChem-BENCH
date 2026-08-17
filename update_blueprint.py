import re
import sys

blueprint_path = r"C:\Users\ansac\.gemini\antigravity\brain\cdbceb39-4e68-416c-a248-6dea9eab8f92\CoChem-BENCH_File_Blueprint.md"

def update_status(filename, status):
    with open(blueprint_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Escape filename for regex
    escaped_filename = re.escape(filename)
    # The file path in the blueprint has backslashes, so we match it exactly
    pattern = r'- \[(?: |/)\] (\[[A-Z]+\]) ' + escaped_filename + r'\b'
    replacement = r'- [' + status + r'] \1 ' + filename.replace('\\', r'\\')

    new_content, count = re.subn(pattern, replacement, content)
    if count == 0:
        print(f"Could not find or update: {filename}")
    else:
        with open(blueprint_path, 'w', encoding='utf-8') as f:
            f.write(new_content)
        print(f"Updated {filename} to [{status}]")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        update_status(sys.argv[1], sys.argv[2])

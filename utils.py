import subprocess
import os

def create_file(filename, content):
    try:
        with open(filename, 'w') as file:
            file.write(content)
        print(f"✅ File '{filename}' created successfully. \n")
    except IOError:
        print(f"❌ Error: Could not create file '{filename}'. \n")


def check_tool_availability(tool_name):
    try:
        subprocess.run([tool_name, "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False
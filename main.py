import sys
import time
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run Generator or AutoTool based on the provided flag.")
    parser.add_argument("file_path", help="The path to the input file.")
    group = parser.add_mutually_exclusive_group(required=True) #just one, either gen or tooq
    group.add_argument("-generator", action="store_true", help="Run the generator script.")
    group.add_argument("-gentool", action="store_true", help="Run the gentool script.")
    parser.add_argument("-o", "--output", help="The name of the output file.")
    parser.add_argument("-debug", action="store_true", help="Enable debug mode.", default=False)

    args = parser.parse_args()

    file_path = args.file_path
    output_file = args.output
    debug = args.debug

    if args.generator:
        print("\nRunning generation process!")
        command = ["python3", "generator.py", file_path]
        if output_file:
            command += ["-o", output_file]
        if debug:
            command.append("-debug")

    elif args.gentool:
        print("\nRunning tool process!")
        command = ["python3", "gen_tool.py", file_path]

    else:
        print("You must specify either -generator or -gentool.")
        sys.exit(1)

    try:
        start_time = time.time()
        subprocess.run(command, check=True)
        end_time = time.time()
        elapsed_time = round(end_time - start_time, 2)
        print(f"\nThe process took {elapsed_time} seconds!")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the script: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

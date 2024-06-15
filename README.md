# Chemical Reaction Generator 🧪

## Overview ℹ️
This tool allows for the generation of chemical species and reactions based on user input. The tool has two main modes: `generator` and `gentool`, which can be selected using the appropriate flags.

## Usage 🚀

### Command Line Arguments 📋

- `file_path`: The path to the input file containing the chemical data (how the file is supposed to look in order to work is shown in test folder, or in documentation relative to the mode you're using in docs folder).
- `-o`, `--output`: The name of the output file where the generated data will be saved.
- `-debug`: Enable debug mode to print additional information (only for **generator mode**).

### Flags 🚩

- `-generator`: Use the generator mode for generating species and reactions.
- `-gentool`: Use the gentool mode as a quick way to generate an input file for the actual generator (check the relative docs).

### Example 🌟

To run the tool in generator mode:
```bash
python main.py test/input/input.txt -generator -o test/output/output.txt -debug
```

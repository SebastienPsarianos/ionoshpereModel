import argparse
import json
import re
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convert Hardy Hall file to pivoted JSON.")
    parser.add_argument('-i', '--input', required=True,
                        help="Path to the input text file")
    parser.add_argument('-o', '--output', required=True,
                        help="Path to the output JSON file")
    return parser.parse_args()


def clean_and_parse(content):
    # 1. Clean [source] tags
    # We use a raw string (r'') here. This tells Python to treat backslashes literally.
    # Pattern explanation: Matches ""
    content = re.sub(r'\\', '', content)

    # 2. Regex to capture the 6 columns
    # We use raw strings (r'') concatenated inside parentheses.
    # This is the safest way to write multi-line regex in Python.
    pattern = re.compile(
        r'(K\d+)\s*,\s*'          # Group 1: Kp (e.g. K0)
        r'([A-Za-z0-9]+)\s*,\s*'  # Group 2: Term (e.g. Cos1)
        r'([-\d.]+)\s*,\s*'       # Group 3: max_value
        r'([-\d.]+)\s*,\s*'       # Group 4: max_latitude
        r'([-\d.]+)\s*,\s*'       # Group 5: up_slope
        r'([-\d.]+)'              # Group 6: down_slope
    )

    temp_data = {}
    properties = ["max_value", "max_latitude", "up_slope", "down_slope"]

    for match in pattern.finditer(content):
        kp = match.group(1)
        term = match.group(2)

        values = {
            "max_value": float(match.group(3)),
            "max_latitude": float(match.group(4)),
            "up_slope": float(match.group(5)),
            "down_slope": float(match.group(6))
        }

        # Initialize Kp structure if missing
        if kp not in temp_data:
            temp_data[kp] = {}
            for prop in properties:
                temp_data[kp][prop] = {"const": None, "cos": [], "sin": []}

        # Route data to the correct bucket
        if term == "Constant":
            for prop in properties:
                temp_data[kp][prop]["const"] = values[prop]

        elif term.startswith("Cos"):
            idx = int(term.replace("Cos", ""))
            for prop in properties:
                temp_data[kp][prop]["cos"].append((idx, values[prop]))

        elif term.startswith("Sin"):
            idx = int(term.replace("Sin", ""))
            for prop in properties:
                temp_data[kp][prop]["sin"].append((idx, values[prop]))

    # 3. Finalize: Sort lists and clean up
    final_output = {}
    for kp in sorted(temp_data.keys()):
        final_output[kp] = {}
        for prop in properties:
            final_output[kp][prop] = {}

            # Constant
            final_output[kp][prop]["const"] = temp_data[kp][prop]["const"]

            # Cosine (sort by index, then store only the value)
            sorted_cos = sorted(temp_data[kp][prop]["cos"], key=lambda x: x[0])
            final_output[kp][prop]["cos"] = [val for idx, val in sorted_cos]

            # Sine (sort by index, then store only the value)
            sorted_sin = sorted(temp_data[kp][prop]["sin"], key=lambda x: x[0])
            final_output[kp][prop]["sin"] = [val for idx, val in sorted_sin]

    return final_output


def main():
    args = parse_arguments()

    try:
        with open(args.input, 'r', encoding='utf-8') as f:
            raw_content = f.read()
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)

    structured_data = clean_and_parse(raw_content)

    try:
        with open(args.output, 'w', encoding='utf-8') as f:
            json.dump(structured_data, f, indent=4)
        print(f"Successfully processed {len(structured_data)} Kp groups.")
        print(f"Output saved to: {args.output}")
    except IOError as e:
        print(f"Error writing to file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

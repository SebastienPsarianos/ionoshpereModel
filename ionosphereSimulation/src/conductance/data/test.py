import re

# Input/Output filenames
INPUT_FILE = "hardy_hall_coefficients.txt"
OUTPUT_FILE = "HardyHallCoefficients.cpp"


def parse_hardy_file(filepath):
    # Data structure: data[kp_index][term_type (0-3)] = { 'const': val, 'cos': [6], 'sin': [6] }
    data = {}

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f):
            parts = [p.strip() for p in line.split(',')]

            # Skip invalid lines or headers
            if not parts or len(parts) < 6:
                continue
            if parts[0] == "Kp" or "term" in parts[1]:
                continue
            if not parts[0].startswith('K'):
                continue

            try:
                # Extract Kp index (e.g., "K0" -> 0)
                kp = int(parts[0][1])
            except ValueError:
                continue

            if kp not in data:
                data[kp] = [{'c': 0.0, 'cos': [0.0]*6, 'sin': [0.0]*6}
                            for _ in range(4)]

            label = parts[1]
            try:
                values = [float(parts[2]), float(parts[3]),
                          float(parts[4]), float(parts[5])]
            except ValueError:
                continue

            # Map coefficients to the correct array
            if label == "Constant":
                for i in range(4):
                    data[kp][i]['c'] = values[i]
            elif label.startswith("Cos"):
                try:
                    n = int(label[3:]) - 1
                    if 0 <= n < 6:
                        for i in range(4):
                            data[kp][i]['cos'][n] = values[i]
                except ValueError:
                    pass
            elif label.startswith("Sin"):
                try:
                    n = int(label[3:]) - 1
                    if 0 <= n < 6:
                        for i in range(4):
                            data[kp][i]['sin'][n] = values[i]
                except ValueError:
                    pass

    return data


def write_cpp(data):
    with open(OUTPUT_FILE, 'w') as f:
        f.write('#include "HardyConductance.hpp"\n\n')
        f.write('namespace ionosphere {\n\n')
        f.write('// Auto-generated coefficients from Hardy et al. 1987\n')

        # Outer array (std::array) needs double braces {{ ... }}
        f.write('const std::array<EpsteinParams, 7> HARDY_COEFFS = {{\n')

        for kp in sorted(data.keys()):
            f.write(f'    // Kp = {kp}\n')

            # INNER STRUCT (EpsteinParams): Needs SINGLE brace { ... }
            f.write('    {\n')

            for param in data[kp]:
                f.write('        { ')  # Start FourierCoeffs struct
                f.write(f'{param["c"]:<10.7f}, ')

                # std::array inside struct: needs braces
                f.write(
                    '{' + ', '.join(f'{x:<10.7f}' for x in param['cos']) + '}, ')
                f.write(
                    '{' + ', '.join(f'{x:<10.7f}' for x in param['sin']) + '}')

                f.write(' },\n')  # End FourierCoeffs

            f.write('    },\n')  # End EpsteinParams (Single brace!)

        f.write('}};\n\n')  # End outer std::array (Double braces)
        f.write('} // namespace ionosphere\n')


if __name__ == "__main__":
    try:
        data = parse_hardy_file(INPUT_FILE)
        if not data:
            print("Error: No data found.")
        else:
            write_cpp(data)
            print(f"Successfully generated {OUTPUT_FILE}")
    except Exception as e:
        print(f"Error: {e}")

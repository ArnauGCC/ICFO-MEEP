#!/usr/bin/env python3
"""
import re
import argparse


def extract_matrix(filename, parameter):

    # Data structure:
    # {
    #     height1: {index1: value, index2: value, ...},
    #     height2: {index1: value, index2: value, ...},
    # }
    data = {}

    current_height = None
    current_index = None

    # Matches:
    # Height: 50.0nm
    height_re = re.compile(r"^\s*Height:\s*([-+0-9.eE]+)\s*nm")

    # Matches:
    # INDEX: 1.05:
    index_re = re.compile(r"^\s*INDEX:\s*([-+0-9.eE]+)\s*:")

    # Matches:
    # PR: 0.0051104578983831314
    parameter_re = re.compile(
        rf"^\s*{re.escape(parameter)}\s*:\s*([-+0-9.eE]+)"
    )

    with open(filename, "r") as f:
        for line in f:

            # New height
            match = height_re.match(line)
            if match:
                current_height = float(match.group(1))
                data[current_height] = {}
                current_index = None
                continue

            # New index
            match = index_re.match(line)
            if match:
                if current_height is None:
                    continue

                current_index = float(match.group(1))
                continue

            # Requested parameter
            match = parameter_re.match(line)
            if match:
                if current_height is None or current_index is None:
                    continue

                value = float(match.group(1))
                data[current_height][current_index] = value

    if not data:
        raise ValueError("No Height blocks found in the input file.")

    # Collect all indices and sort them
    indices = sorted({
        index
        for height_data in data.values()
        for index in height_data
    })

    heights = sorted(data.keys())

    # Build matrix
    matrix = []

    for height in heights:
        row = []

        for index in indices:
            if index not in data[height]:
                row.append(None)
            else:
                row.append(data[height][index])

        matrix.append(row)

    return matrix, heights, indices


def main():
    parser = argparse.ArgumentParser(
        description="Extract a parameter matrix from simulation output."
    )

    parser.add_argument(
        "input_file",
        help="Input simulation output file"
    )

    parser.add_argument(
        "parameter",
        help="Parameter to extract, e.g. PR, PL, 0L, 1R"
    )

    args = parser.parse_args()

    matrix, heights, indices = extract_matrix(
        args.input_file,
        args.parameter
    )

    print('[')
    for row in matrix:
        print(row, ',')
    print(']')


if __name__ == "__main__":
    main()
"""

import re
import sys


def extract_parameter(filename, parameter):
    """
    Extract all values corresponding to `parameter` from the file.

    Example:
        extract_parameter("results.txt", "PR")
    """
    pattern = rf"^\s*{re.escape(parameter)}:\s*([-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?)"

    values = []

    with open(filename, "r") as f:
        for line in f:
            match = re.match(pattern, line)
            if match:
                values.append(float(match.group(1)))

    return values


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_parameter.py <input_file> <parameter>")
        sys.exit(1)

    filename = sys.argv[1]
    parameter = sys.argv[2]

    values = extract_parameter(filename, parameter)

    print(f"{parameter} = {values}")
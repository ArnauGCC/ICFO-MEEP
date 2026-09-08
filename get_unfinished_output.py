#!/usr/bin/env python3

import re
import sys
import numpy as np


def extract_matrix(filename, parameter):
    """
    Extract a parameter from a simulation output file.

    Rows    = grat_depth_factor
    Columns = grat_duty_cycle
    """

    depth_blocks = {}
    current_depth = None
    current_duty = None

    with open(filename, "r") as f:

        for line in f:

            # ---------------------------------------------------------
            # Match grat_depth_factor
            #
            # Example:
            # Executing with grat_depth_factor:   0.025
            # ---------------------------------------------------------
            match = re.search(
                r"Executing with grat_height:\s*([-+0-9.eE]+)",
                line
            )

            if match:
                current_depth = float(match.group(1))

                if current_depth not in depth_blocks:
                    depth_blocks[current_depth] = {}

                current_duty = None
                continue

            # ---------------------------------------------------------
            # Match grat_duty_cycle
            #
            # Example:
            # Executing with grat_duty_cycle:   0.025
            # ---------------------------------------------------------
            match = re.search(
                r"Executing with grat_depth_factor:\s*([-+0-9.eE]+)",
                line
            )

            if match:
                current_duty = float(match.group(1))
                continue

            # ---------------------------------------------------------
            # Match requested parameter
            #
            # Examples:
            # PR: 0.0010925142429455018
            # PL: 0.00014912595433007615
            # 0L: 1.0737228423177608e-12
            # ---------------------------------------------------------
            match = re.match(
                rf"\s*{re.escape(parameter)}:\s*([-+0-9.eE]+)",
                line
            )

            if match and current_depth is not None and current_duty is not None:

                value = float(match.group(1))

                depth_blocks[current_depth][current_duty] = value

    # -------------------------------------------------------------
    # Sort depth factors -> matrix rows
    # Sort duty cycles   -> matrix columns
    # -------------------------------------------------------------
    depths = sorted(depth_blocks.keys())

    duties = sorted({
        duty
        for block in depth_blocks.values()
        for duty in block.keys()
    })

    # -------------------------------------------------------------
    # Create matrix
    # -------------------------------------------------------------
    matrix = np.full(
        (len(depths), len(duties)),
        np.nan,
        dtype=float
    )

    # -------------------------------------------------------------
    # Fill matrix
    # -------------------------------------------------------------
    for i, depth in enumerate(depths):

        for j, duty in enumerate(duties):

            if duty in depth_blocks[depth]:
                matrix[i, j] = depth_blocks[depth][duty]

    return depths, duties, matrix


def print_matrix(matrix):
    """
    Print matrix with:
        - one row per line
        - comma after every value
        - comma after every row

    Example:

    [
        [1.0, 2.0, 3.0, ],
        [4.0, 5.0, 6.0, ],
    ]
    """

    print("[")

    for row in matrix:

        print("    [", end="")

        for value in row:
            print(f"{value}, ", end="")

        print("],")

    print("]")


# =================================================================
# Main
# =================================================================

if __name__ == "__main__":

    # -------------------------------------------------------------
    # Check command line arguments
    # -------------------------------------------------------------
    if len(sys.argv) != 3:

        print(
            f"Usage: {sys.argv[0]} <input_file> <parameter>"
        )

        print(
            f"Example: {sys.argv[0]} output.txt PR"
        )

        sys.exit(1)

    # -------------------------------------------------------------
    # Get arguments
    # -------------------------------------------------------------
    filename = sys.argv[1]
    parameter = sys.argv[2]

    # -------------------------------------------------------------
    # Extract matrix
    # -------------------------------------------------------------
    depths, duties, matrix = extract_matrix(
        filename,
        parameter
    )

    # -------------------------------------------------------------
    # Print information
    # -------------------------------------------------------------
    print(f"\nParameter: {parameter}")
    print(f"Number of depth factors: {len(depths)}")
    print(f"Number of duty cycles:   {len(duties)}")

    print("\nMatrix:")

    # -------------------------------------------------------------
    # Print matrix
    # -------------------------------------------------------------
    print_matrix(matrix)

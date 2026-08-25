import re


def extract_quantity(text, quantity):
    """
    Extract `quantity` (e.g. 0L, 0R, PR, PL, 1L) from every
    WG WIDTH cycle.

    Returns:
        list[list[float]]

    Rows    = WG WIDTH cycles
    Columns = frequency points
    Missing values are filled with 0.0
    """

    # ------------------------------------------------------------
    # Find the first WG WIDTH
    # Everything before it is ignored.
    # ------------------------------------------------------------
    first_cycle = re.search(
        r"WG WIDTH:\s*[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?",
        text
    )

    if not first_cycle:
        return []

    text = text[first_cycle.start():]

    # ------------------------------------------------------------
    # Split into WG WIDTH cycles
    # ------------------------------------------------------------
    cycles = re.split(
        r"(?=WG WIDTH:\s*[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
        text
    )

    cycles = [
        cycle.strip()
        for cycle in cycles
        if cycle.strip()
    ]

    results = []

    # ------------------------------------------------------------
    # Process every WG WIDTH cycle
    # ------------------------------------------------------------
    for cycle in cycles:

        # Find every frequency block
        freq_pattern = (
            r"Freq:\s*"
            r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
            r"\s*:\s*"
            r"(.*?)(?=\n\s*Freq:|\Z)"
        )

        freq_matches = re.findall(
            freq_pattern,
            cycle,
            flags=re.DOTALL
        )

        row = []

        for freq, block in freq_matches:

            # Match:
            #
            # PR: 0.123
            #
            # or:
            #
            # 1L: np.float64(0.123)
            #
            value_pattern = (
                rf"^\s*{re.escape(quantity)}\s*:\s*"
                rf"(?:np\.float64\()? "
                rf"([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)"
            )

            match = re.search(
                value_pattern,
                block,
                flags=re.MULTILINE | re.VERBOSE
            )

            if match:
                row.append(float(match.group(1)))
            else:
                row.append(0.0)

        results.append(row)

    # ------------------------------------------------------------
    # Make every row have the same number of columns
    # ------------------------------------------------------------
    max_columns = max(
        (len(row) for row in results),
        default=0
    )

    for row in results:
        while len(row) < max_columns:
            row.append(0.0)

    return results


# ================================================================
# Example
# ================================================================

if __name__ == "__main__":

    with open("outs.txt", "r", encoding="utf-8") as f:
        text = f.read()

    quantity = input(
        "Quantity to extract (e.g. 0L, 0R, PR, PL, 1L): "
    ).strip()

    matrix = extract_quantity(text, quantity)

    # Print each row separately with a line return
    for row in matrix:
        print(row, ',')
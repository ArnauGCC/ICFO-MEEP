import re


def extract_data(filename, parameter):
    """
    Extract a parameter from a simulation output file.

    Returns:
        A Python list of lists.
        Each row corresponds to one WG WIDTH.
        Each column corresponds to one frequency.
    """

    parameter = str(parameter).strip().upper()

    #if parameter not in {"0", "1", "2", "P", "S"}:
    #    raise ValueError("Parameter must be one of: 0, 1, 2, P, S")

    with open(filename, "r", encoding="utf-8") as file:
        text = file.read()

    # Split the file into WG WIDTH blocks
    blocks = re.split(r"(?=WG WIDTH:\s*)", text)

    matrix = []

    for block in blocks:

        # Find WG WIDTH
        width_match = re.search(
            r"WG WIDTH:\s*([0-9.eE+-]+)",
            block
        )

        if width_match is None:
            continue

        values = []

        # Find each Freq section
        freq_blocks = re.finditer(
            r"Freq:\s*[0-9.eE+-]+:\s*(.*?)(?=Freq:|$)",
            block,
            re.DOTALL
        )

        for match in freq_blocks:

            freq_data = match.group(1)

            # Find the requested parameter
            value_match = re.search(
                rf"^\s*{re.escape(parameter)}:\s*"
                r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
                freq_data,
                re.MULTILINE
            )

            if value_match:
                values.append(float(value_match.group(1)))

        if values:
            matrix.append(values)

    if not matrix:
        return []

    # Check that all rows have the same length
    lengths = [len(row) for row in matrix]

    if len(set(lengths)) != 1:
        raise ValueError(
            f"Different numbers of frequencies found: {lengths}"
        )

    return matrix


def main():

    filename = "server/output.txt".strip()

    parameter = '0L'

    try:

        matrix = extract_data(filename, parameter)

        print("\nMatrix:")
        print("[")

        for i, row in enumerate(matrix):

            if i < len(matrix) - 1:
                print(f"    {row},")
            else:
                print(f"    {row}")

        print("]")

    except FileNotFoundError:
        print(f"File not found: {filename}")

    except ValueError as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
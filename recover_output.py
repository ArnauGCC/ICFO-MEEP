import numpy as np
import ast

with open("server/output.txt") as f:
    lines = [line.rstrip() for line in f]

freqs = []
eff = []
parameters = []

current_steps = None
current_eff = None

i = 0
while i < len(lines):

    line = lines[i]

    # A new frequency starts -> save the previous one
    if line.startswith("FREQUENCY"):
        if current_steps is not None:
            idx = np.argmax(current_eff)
            eff.append(current_eff[idx])
            parameters.append(current_steps[idx])

        freqs.append(float(line.split("=")[1]))
        current_steps = None
        current_eff = None

    elif line == "STEPS:":
        i += 1
        text = lines[i]
        while "]" not in text:
            i += 1
            text += " " + lines[i]

        current_steps = np.fromstring(
            text.replace("[", "").replace("]", ""),
            sep=" "
        )

    elif line == "EFF:":
        i += 1
        text = lines[i]
        while "]" not in text:
            i += 1
            text += " " + lines[i]

        current_eff = np.array(ast.literal_eval(text))

    i += 1

# Save the last frequency
if current_steps is not None:
    idx = np.argmax(current_eff)
    eff.append(current_eff[idx])
    parameters.append(current_steps[idx])

print(freqs)
print(eff)
print(parameters)
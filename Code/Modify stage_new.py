import re
import numpy as np

INPUT_FILE = "Other files/stage_new.dat"
OUTPUT_FILE = "Output/stage_new modified.dat"

# number of points in block2/block3 for each blade row
NPTS = {
    1: 111,
    3: 106
}


def numbers_from_lines(lines):
    vals = []
    for line in lines:
        vals.extend(float(x) for x in line.split())
    return np.array(vals)


def format_block(arr):
    out = []
    for i in range(0, len(arr), 8):
        row = arr[i:i+8]
        out.append(" " + " ".join(f"{x: .6f}" for x in row))
    return out


with open(INPUT_FILE) as f:
    lines = f.readlines()

i = 0
current_row = None

# ---------- GLOBAL CHANGES ----------
for idx, line in enumerate(lines):

    # NSTEPS_MAX
    if "NSTEPS_MAX, CONLIM" in line:
        lines[idx+1] = "     20000  0.000100\n"

    # RFMIX/FEXTRAP line
    if "RFMIX,    FEXTRAP,   FSMTHB,    FANGLE" in line:
        lines[idx+1] = "  0.025000  0.800000  1.000000  0.800000\n"

    # IPOUT SFEXIT NSFEXIT
    if "IPOUT  SFEXIT  NSFEXIT" in line:
        lines[idx+1] = "    1  0.100000    10\n"

    # ILOS NLOS IBOUND
    if "ILOS" in line and "NLOS" in line and "IBOUND" in line:
        lines[idx+1] = "        100        5         0\n"

    # Marker block
    if "MARKER FOR VARIABLES TO BE SENT TO THE OUTPUT FILE." in line:
        lines[idx+1] = " 2 2 2 2 2 2 2 2 0 0 2 0 0\n"

    # Stream surfaces block
    if "STREAM SURFACES ON WHICH RESULTS ARE TO BE SENT TO" in line:
        lines[idx+1] = (
            " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
            "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"
        )

# Blade row changes
while i < len(lines):

    line = lines[i]

    m = re.search(r'BLADE ROW NUMBER\s*=\s*(\d+)', line)
    if m:
        current_row = int(m.group(1))
        print("Entered blade row", current_row)

        # ---------- TIP CLEARANCE MODIFICATIONS ----------
        if current_row == 2:

            for k in range(i, min(i+50, len(lines))):
                if "KTIPSTART" in lines[k]:

                    lines[k+1] = "         34        37\n"

                    extra = [
                        "      FRACTIP1  FRACTIP2\n",
                        "       0.0103    0.0103\n",
                        "      FTHICK(K) K=1,KM\n",
                        "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 "
                        "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 "
                        "1.0 0.9 0.5 0.0 0.0 0.0 0.0\n"
                    ]

                    lines[k+2:k+2] = extra
                    break


        if current_row == 4:

            for k in range(i, min(i+50, len(lines))):
                if "KTIPSTART" in lines[k]:

                    lines[k+1] = "         35        37\n"

                    extra = [
                        "      FRACTIP1  FRACTIP2\n",
                        "       0.0059    0.0059\n",
                        "      FTHICK(K) K=1,KM\n",
                        "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 "
                        "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 "
                        "1.0 1.0 0.9 0.5 0.0 0.0 0.0\n"
                    ]

                    lines[k+2:k+2] = extra
                    break

    # Addition
    if current_row in [1, 3]:

        # start of section?
        vals = line.split()
        if vals == ['1.00000', '0.00000', '0']:
            print(f"Found section in blade row {current_row}")

            npts = NPTS[current_row]

            # ----- skip section header and block1 -----
            j = i + 1

            count = 0
            while count < npts:
                count += len(lines[j].split())
                j += 1

            # now j points to line: "1.000000 0.000000"
            sep1 = j
            j += 1

            # ----- read block2 -----
            block2_start = j
            count = 0
            while count < npts:
                count += len(lines[j].split())
                j += 1
            block2_end = j

            # separator line containing only 1.000000
            sep2 = j
            j += 1

            # ----- read block3 -----
            block3_start = j
            count = 0
            while count < npts:
                count += len(lines[j].split())
                j += 1
            block3_end = j

            # compute sum
            block2 = numbers_from_lines(lines[block2_start:block2_end])
            block3 = numbers_from_lines(lines[block3_start:block3_end])

            summed = block2 + block3
            new_block2 = format_block(summed)

            # replace block2 only
            lines[block2_start:block2_end] = [s + "\n" for s in new_block2]

            # adjust index after replacement
            i = block3_end
            continue

    i += 1


with open(OUTPUT_FILE, "w") as f:
    f.writelines(lines)

print("Written to", OUTPUT_FILE)
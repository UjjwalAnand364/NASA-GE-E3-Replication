import re
import numpy as np

INPUT_FILE = "stage_new.dat"
OUTPUT_FILE = "C:\cygwin64\Multall_Codes\stage_new.dat"

NPTS = {1:111, 3:106}

ANGLE_BLOCKS = {
    1:[
        "N_ANGLES\n","3\n",
        "0.0  0.5  1.0\n",
        "-21.0  -15.0  -12.5\n",
        "67.5  70.0  67.0\n",
        "67.5  70.0  67.0\n"
    ],
    2:[
        "N_ANGLES\n","3\n",
        "0.0  0.5  1.0\n",
        "30.0  32.5  32.0\n",
        "-59.0  -59.0  -58.0\n",
        "-59.0  -59.0  -58.0\n"
    ],
    3:[
        "N_ANGLES\n","3\n",
        "0.0  0.5  1.0\n",
        "-14.0  -16.0  -15.0\n",
        "64.0  59.0  64.5\n",
        "64.0  59.0  64.5\n"
    ],
    4:[
        "N_ANGLES\n","3\n",
        "0.0  0.5  1.0\n",
        "32.0  13.0  12.0\n",
        "-52.0  -54.5  -55.0\n",
        "-52.0  -54.5  -55.0\n"
    ]
}



def numbers_from_lines(lines):
    vals=[]
    for line in lines:
        vals.extend(map(float,line.split()))
    return np.array(vals)


def format_block(arr):
    out=[]
    for i in range(0,len(arr),8):
        out.append(" "+" ".join(f"{x: .6f}" for x in arr[i:i+8]))
    return out


with open(INPUT_FILE) as f:
    lines=f.readlines()


# ==========================================================
# GLOBAL CHANGES
# ==========================================================

for idx,line in enumerate(lines):

    if "NSTEPS_MAX, CONLIM" in line:
        lines[idx+1]="     20000  0.000100\n"

    elif "RFMIX,    FEXTRAP,   FSMTHB,    FANGLE" in line:
        lines[idx+1]="  0.010000  0.000000  1.000000  0.800000\n"

    elif "IPOUT  SFEXIT  NSFEXIT" in line:
        lines[idx+1]="    1  0.100000    10\n"

    elif "ILOS" in line and "NLOS" in line:
        lines[idx+1]="        100        5         0\n"

    elif "MARKER FOR VARIABLES TO BE SENT TO THE OUTPUT FILE." in line:
        lines[idx+1]=" 2 2 2 2 2 2 2 2 0 0 2 0 0\n"

    elif "STREAM SURFACES ON WHICH RESULTS ARE TO BE SENT TO" in line:
        lines[idx+1]=" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n"

    elif "MIXING LENGTH LIMITS ON ALL BLADE ROWS" in line:
        new_block = "  0.020000  0.030000  0.040000  0.050000  0.000000  0.500000\n"
        lines[idx+1:idx+5] = [new_block] * 4


# ==========================================================
# PROCESS BLADE ROWS
# ==========================================================

current_row = None
i = 0

while i < len(lines):

    # ----------------------------------------------------------
    # Detect start of a blade row
    # ----------------------------------------------------------
    m = re.search(r'BLADE ROW NUMBER\s*=\s*(\d+)', lines[i])

    if m:
        current_row = int(m.group(1))
        print("Entered blade row", current_row)

        # --------------------------------------------------
        # Search within current blade row header
        # --------------------------------------------------

        k = i

        while k < len(lines):

            if (
                k != i
                and "************STARTING THE INPUT FOR EACH BLADE ROW" in lines[k]
            ):
                break

            # ---------------- IFANGLES ----------------

            if "IF_CUSP_OUT" in lines[k]:
                vals = lines[k+1].split()
                vals[1] = "1"
                lines[k+1] = f"{int(vals[0]):10d}{int(vals[1]):10d}\n"

            # ---------------- Row 2 tip clearance ----------------

            if current_row == 2 and "KTIPSTART" in lines[k]:

                lines[k+1] = "         34        37\n"

                lines[k+2:k+2] = [
                    "      FRACTIP1  FRACTIP2\n",
                    "       0.0103    0.0103\n",
                    "      FTHICK(K) K=1,KM\n",
                    "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 "
                    "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 "
                    "1.0 0.9 0.5 0.0 0.0 0.0 0.0\n"
                ]

                k += 4

            # ---------------- Row 4 tip clearance ----------------

            if current_row == 4 and "KTIPSTART" in lines[k]:

                lines[k+1] = "         35        37\n"

                lines[k+2:k+2] = [
                    "      FRACTIP1  FRACTIP2\n",
                    "       0.0059    0.0059\n",
                    "      FTHICK(K) K=1,KM\n",
                    "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 "
                    "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 "
                    "1.0 1.0 0.9 0.5 0.0 0.0 0.0\n"
                ]

                k += 4

            k += 1

    # ----------------------------------------------------------
    # block2 += block3
    # ----------------------------------------------------------

    if current_row in (1, 3):

        if lines[i].split() == ['1.00000', '0.00000', '0']:

            n = NPTS[current_row]

            j = i + 1

            c = 0
            while c < n:
                c += len(lines[j].split())
                j += 1

            j += 1

            b2s = j

            c = 0
            while c < n:
                c += len(lines[j].split())
                j += 1

            b2e = j

            j += 1

            b3s = j

            c = 0
            while c < n:
                c += len(lines[j].split())
                j += 1

            b3e = j

            b2 = numbers_from_lines(lines[b2s:b2e])
            b3 = numbers_from_lines(lines[b3s:b3e])

            lines[b2s:b2e] = [x + "\n" for x in format_block(b2 + b3)]

            i = b3e
            continue

    i += 1

# ==========================================================
# Insert N_ANGLES blocks
# ==========================================================

new_lines = []
row = 1

i = 0
while i < len(lines):

    # Detect the pair:
    # ***************************************************************
    # ************STARTING THE INPUT FOR EACH BLADE ROW**************
    if (
        i + 1 < len(lines)
        and lines[i].strip() == "***************************************************************"
        and "STARTING THE INPUT FOR EACH BLADE ROW" in lines[i + 1]
    ):

        # Insert previous row's angles before the separator
        if row > 1:
            new_lines.extend(ANGLE_BLOCKS[row - 1])

        new_lines.append(lines[i])
        new_lines.append(lines[i + 1])

        row += 1
        i += 2
        continue

    # Final blade row
    if "STARTING INLET BOUNDARY CONDITION DATA" in lines[i]:
        new_lines.extend(ANGLE_BLOCKS[4])

    new_lines.append(lines[i])
    i += 1

lines = new_lines


with open(OUTPUT_FILE,"w") as f:
    f.writelines(lines)

print("Done.")
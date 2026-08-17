import numpy as np
import sys

def load_blade(filepath):
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]

def find_le_te(x, y):
    le = np.array([(x[0] + x[-1]) / 2, (y[0] + y[-1]) / 2])
    i_te = np.argmax(np.hypot(x - le[0], y - le[1]))
    te = np.array([x[i_te], y[i_te]])
    return le, te, i_te

def convert_to_stagen(dat_file, rotn=0.0, nxpts=200, ifclock=0, ifrev=0,
                       xrot=0.5, yrot=0.5, xcup=0.5, xcdwn=0.25):

    x, y = load_blade(dat_file)
    le, te, i_te = find_le_te(x, y)
    chord_len = np.linalg.norm(te - le)

    # Translate LE to origin, keep true machine orientation
    xc = x - le[0]
    yc = y - le[1]

    xc = xc * 2.54
    yc = yc * 2.54

    # LE and TE metal angles from raw surface tangents (averaged SS+PS)
    n_avg = 5
    le_tan_ss = np.array([x[n_avg]    - x[0],    y[n_avg]    - y[0]])
    le_tan_ps = np.array([x[-n_avg]   - x[-1],   y[-n_avg]   - y[-1]])
    le_tan_ss /= np.linalg.norm(le_tan_ss)
    le_tan_ps /= np.linalg.norm(le_tan_ps)
    le_tan = le_tan_ss + le_tan_ps
    le_tan /= np.linalg.norm(le_tan)

    te_tan_ss = np.array([x[i_te] - x[i_te - n_avg], y[i_te] - y[i_te - n_avg]])
    te_tan_ps = np.array([x[i_te] - x[i_te + n_avg], y[i_te] - y[i_te + n_avg]])
    te_tan_ss /= np.linalg.norm(te_tan_ss)
    te_tan_ps /= np.linalg.norm(te_tan_ps)
    te_tan = te_tan_ss + te_tan_ps
    te_tan /= np.linalg.norm(te_tan)

    betup = np.degrees(np.arctan2(le_tan[1], le_tan[0]))
    betdn = np.degrees(np.arctan2(te_tan[1], te_tan[0]))

    npoints = len(xc)
    lines = []
    lines.append(f"    0                    INTYPE- TYPE OF BLADE GEOMETRY INPUT")
    lines.append(f"  {npoints}  {nxpts}  {ifclock}  {ifrev}          NPOINTS, NXPTS, IFCLOCK, IFREV")
    for bx, by in zip(xc, yc):
        lines.append(f"  {bx:.6f}  {by:.6f}")
    lines.append(f"  0.0000  {xrot:.4f}  {yrot:.4f}                    ROTN,XROT,YROT")
    lines.append(f"    {xcup:.4f}    {xcdwn:.4f}    {betup:.4f}    {betdn:.4f}          XCUP, XCDWN, BETUP, BETDWN")

    return "\n".join(lines), chord_len, le, te, betup, betdn


if __name__ == '__main__':
    import os
    dat_file = sys.argv[1] if len(sys.argv) > 1 else 'r1s1.dat'

    block, chord_len, le, te, betup, betdn = convert_to_stagen(dat_file)

    print(f"chord_len={chord_len:.4f}  LE={le}  TE={te}")
    print(f"BETUP={betup:.4f}  BETDN={betdn:.4f}")

    base = os.path.splitext(os.path.basename(dat_file))[0]
    outfile = 'Final dat files/' + base + '_norm.dat'
    with open(outfile, 'w', newline='') as f:
        f.write(block)
    print(f"Written: {outfile}  ({block.count(chr(10))+1} lines)")
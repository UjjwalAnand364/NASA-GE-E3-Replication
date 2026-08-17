import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import sys

def load_blade(filepath):
    data = np.loadtxt(filepath)
    return data[:, 0], data[:, 1]

def find_le_te(x, y):
    # Data is a closed loop: start and end points are both at the TE
    # (they nearly coincide). LE = point farthest from TE.
    te = np.array([(x[0] + x[-1]) / 2, (y[0] + y[-1]) / 2])
    i_te = np.argmax(np.hypot(x - te[0], y - te[1]))
    te = np.array([x[i_te], y[i_te]])
    le = np.array([(x[0] + x[-1]) / 2, (y[0] + y[-1]) / 2])
    return le, te, i_te

def chord_frame(px, py, le, chord_dir, normal_dir):
    pts = np.stack([px - le[0], py - le[1]], axis=1)
    return pts @ chord_dir, pts @ normal_dir

def clean_surface(xc, yc, chord_len):
    mask = (xc >= 0) & (xc <= chord_len)
    xc, yc = xc[mask], yc[mask]
    order = np.argsort(xc)
    return xc[order], yc[order]

def compute_thickness(x, y, n_query=2000):
    le, te, i_te = find_le_te(x, y)

    chord_vec = te - le
    chord_len = np.linalg.norm(chord_vec)
    chord_dir = chord_vec / chord_len

    # Calculate the geometric twist (stagger) angle
    twist_rad = np.arctan2(chord_vec[1], chord_vec[0])
    twist_deg = np.degrees(twist_rad)

    normal_dir = np.array([-chord_dir[1], chord_dir[0]])

    to_cf = lambda px, py: chord_frame(px, py, le, chord_dir, normal_dir)

    # Split at i_te: two surfaces both going from LE-end toward TE
    sA_xc, sA_yc = to_cf(x[:i_te + 1], y[:i_te + 1])
    sB_xc, sB_yc = to_cf(x[i_te:][::-1], y[i_te:][::-1])

    sA_xc, sA_yc = clean_surface(sA_xc, sA_yc, chord_len)
    sB_xc, sB_yc = clean_surface(sB_xc, sB_yc, chord_len)

    xmin = max(sA_xc.min(), sB_xc.min(), 0.005 * chord_len)
    xmax = min(sA_xc.max(), sB_xc.max(), 0.995 * chord_len)
    xq = np.linspace(xmin, xmax, n_query)

    f_A = interp1d(sA_xc, sA_yc, kind='cubic')
    f_B = interp1d(sB_xc, sB_yc, kind='cubic')
    thickness = np.abs(f_A(xq) - f_B(xq))

    idx_max = np.argmax(thickness)
    t_frac = thickness[idx_max] / chord_len
    x_frac = xq[idx_max] / chord_len

    surfaces = dict(
        ss=(sA_xc, sA_yc) if sA_yc.mean() < sB_yc.mean() else (sB_xc, sB_yc),
        ps=(sB_xc, sB_yc) if sA_yc.mean() < sB_yc.mean() else (sA_xc, sA_yc),
        xq=xq, thickness=thickness, chord_len=chord_len,
        le=le, te=te,
        twist_deg=twist_deg 
    )
    return t_frac, x_frac, surfaces

def plot(x, y, surfaces, t_frac, x_frac, outfile=None):
    chord_len = surfaces['chord_len']
    xq = surfaces['xq']
    ss_xc, ss_yc = surfaces['ss']
    ps_xc, ps_yc = surfaces['ps']
    le = surfaces['le']; te = surfaces['te']

    fig, axes = plt.subplots(3, 1, figsize=(10, 12))

    axes[0].plot(x, y, 'b-', lw=1.2)
    axes[0].plot(*le, 'go', ms=8, label='LE')
    axes[0].plot(*te, 'rs', ms=8, label='TE')
    axes[0].plot([le[0], te[0]], [le[1], te[1]], 'k--', lw=1, label='Chord')
    axes[0].set_aspect('equal'); axes[0].grid(True, alpha=0.3)
    axes[0].set_title('Blade profile (original coordinates)'); axes[0].legend()

    axes[1].plot(ss_xc / chord_len, ss_yc / chord_len, 'b-', lw=1.5, label='SS')
    axes[1].plot(ps_xc / chord_len, ps_yc / chord_len, 'r-', lw=1.5, label='PS')
    axes[1].axvline(x_frac, color='k', ls='--', lw=1, label=f'max-t @ x/c={x_frac:.3f}')
    axes[1].set_aspect('equal'); axes[1].grid(True, alpha=0.3)
    axes[1].set_xlabel('x/c'); axes[1].set_ylabel('y/c (chord-normal)')
    axes[1].set_title('Chord-aligned profile'); axes[1].legend()

    axes[2].plot(xq / chord_len, surfaces['thickness'] / chord_len, 'k-', lw=1.5)
    axes[2].axvline(x_frac, color='r', ls='--', lw=1)
    axes[2].axhline(t_frac, color='r', ls='--', lw=1,
                    label=f't/c = {t_frac:.4f}  @  x/c = {x_frac:.4f}')
    axes[2].set_xlabel('x/c'); axes[2].set_ylabel('t/c')
    axes[2].set_title('Thickness distribution (chord-normal)'); axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    if outfile:
        plt.savefig(outfile, dpi=150)
        print(f"Plot saved to {outfile}")
    else:
        plt.show()

if __name__ == '__main__':
    filepath = sys.argv[1] if len(sys.argv) > 1 else 'r1s1.dat'
    x, y = load_blade(filepath)
    t_frac, x_frac, surfaces = compute_thickness(x, y)

    print(f"Chord length            : {surfaces['chord_len']:.6f}")
    print(f"Max thickness / chord   : {t_frac:.4f}")
    print(f"x/c at max thickness    : {x_frac:.4f}")
    print(f"Geometric Twist Angle   : {surfaces['twist_deg']:.4f} degrees")

    # plot(x, y, surfaces, t_frac, x_frac, outfile='thickness_analysis3.png')

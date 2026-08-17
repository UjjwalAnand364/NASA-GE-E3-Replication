import numpy as np
import struct

def read_fortran_record(f):
    marker = f.read(4)
    if len(marker) < 4:
        return None
    nbytes = struct.unpack('i', marker)[0]
    data = f.read(nbytes)
    f.read(4)
    return data


def read_flow_out(filename='Other files/flow_out'):

    IM, JM, KM = 37, 434, 37

    gamma = 1.333
    R_gas = 287.0
    Cp = 1148.86

    with open(filename, 'rb') as f:

        # Record 1
        rec = read_fortran_record(f)
        NSTEP = struct.unpack('i', rec)[0]
        print(f"NSTEPS = {NSTEP}")

        def read_array():

            rec = read_fortran_record(f)

            arr = np.frombuffer(rec, dtype=np.float32)

            expected = IM * JM * KM

            if len(arr) != expected:
                print(f"ERROR: got {len(arr)} values, expected {expected}")

            return arr.reshape((IM, JM, KM), order='F')

        # Actual MULTALL flow_out variables
        DEN     = read_array()     # RO
        DENVX   = read_array()     # ROVX
        DENVR   = read_array()     # ROVR
        DENVT   = read_array()     # ROVT
        DENE    = read_array()     # ROE
        ROSUB   = read_array()
        VISCVAR = read_array()
        SPARE   = read_array()

    print("\nDensity range:")
    print(DEN.min(), DEN.max())

    # Prevent divide by zero
    DEN_safe = np.where(np.abs(DEN) < 1e-10, 1e-10, DEN)

    # Velocities
    Vx = DENVX / DEN_safe
    Vr = DENVR / DEN_safe
    Vt = DENVT / DEN_safe

    V = np.sqrt(Vx**2 + Vr**2 + Vt**2)

    # Internal energy
    E = DENE / DEN_safe

    T_static = (E - 0.5 * V**2) / (Cp - R_gas)

    P_static = DEN_safe * R_gas * T_static

    T_total = T_static + V**2 / (2 * Cp)

    P_total = np.where(
        T_static > 0,
        P_static * (T_total / T_static) ** (gamma / (gamma - 1)),
        0.0
    )

    a = np.sqrt(np.maximum(gamma * R_gas * T_static, 0))

    Mach = np.where(a > 0, V / a, 0)

    # Mid-span
    K_mid = KM // 2

    print("\n=== MID-SPAN RESULTS ===")

    print(f"{'J':>5} {'Mach':>8} {'Pstat(Pa)':>12} {'Ptot(Pa)':>12} {'Ttot(K)':>10}")

    for j in range(JM):

        mach = np.mean(Mach[:, j, K_mid])
        ps   = np.mean(P_static[:, j, K_mid])
        pt   = np.mean(P_total[:, j, K_mid])
        tt   = np.mean(T_total[:, j, K_mid])

        print(f"{j+1:5d} {mach:8.4f} {ps:12.1f} {pt:12.1f} {tt:10.2f}")


if __name__ == '__main__':

    import sys

    fname = sys.argv[1] if len(sys.argv) > 1 else 'Other files/flow_out'

    read_flow_out(fname)
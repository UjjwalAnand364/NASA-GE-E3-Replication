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

def read_flow_out(filename='flow_out'):
    # Grid dimensions - Fortran stores as (K, I, J) column major
    # so reshape as (KM, IM, JM) then transpose to (IM, JM, KM)
    IM, JM, KM = 37, 217, 37
    gamma = 1.333
    R_gas = 287.0
    Cp    = 1148.86

    with open(filename, 'rb') as f:
        # Record 1: NSTEPS
        rec = read_fortran_record(f)
        NSTEPS = struct.unpack('i', rec)[0]
        print(f"NSTEPS = {NSTEPS}")

        def read_array():
            rec = read_fortran_record(f)

            # Use double precision
            arr = np.frombuffer(rec, dtype=np.float64)

            expected = IM * JM * KM
            if len(arr) != expected:
                print(f"Expected {expected} values but got {len(arr)}")
                print(f"Record length = {len(rec)} bytes")

            # MULTALL Fortran ordering
            return arr.reshape((IM, JM, KM), order='F')

        X      = read_array()
        R      = read_array()
        RT     = read_array()
        DEN    = read_array()
        DENVX  = read_array()
        DENVR  = read_array()
        DENRVT = read_array()
        DENE   = read_array()

    # Sanity check on R - should be ~0.3-0.4m for HPT
    print(f"\nR range: {R.min():.4f} to {R.max():.4f} m")
    print(f"X range: {X.min():.4f} to {X.max():.4f} m")
    print(f"DEN range: {DEN.min():.4f} to {DEN.max():.4f} kg/m3")

    # Derived quantities
    # Avoid division by zero
    DEN_safe = np.where(np.abs(DEN) < 1e-10, 1e-10, DEN)
    R_safe   = np.where(np.abs(R)   < 1e-10, 1e-10, R)

    Vx  = DENVX  / DEN_safe
    Vr  = DENVR  / DEN_safe
    Vt  = DENRVT / (DEN_safe * R_safe)
    Vm  = np.sqrt(Vx**2 + Vr**2)
    V   = np.sqrt(Vx**2 + Vr**2 + Vt**2)
    E        = DENE / DEN_safe
    T_static = (E - 0.5*V**2) / (Cp - R_gas)
    P_static = DEN_safe * R_gas * T_static
    T_total  = T_static + V**2 / (2*Cp)
    P_total  = P_static * np.where(T_static > 0,
                (T_total/np.where(T_static>0, T_static, 1))**(gamma/(gamma-1)), 0)
    a        = np.sqrt(np.where(gamma * R_gas * T_static > 0,
                gamma * R_gas * T_static, 0))
    Mach_abs = np.where(a > 0, V/a, 0)
    RVt      = R * Vt

    # Mid-span
    K_mid = KM // 2
    print(f"\n=== FLOW QUANTITIES AT MID-SPAN (K={K_mid+1}) - Pitchwise average ===")
    print(f"{'J':>5} {'Vx':>8} {'Vm':>8} {'Vt':>8} {'Mach':>7} "
          f"{'P_stat(Pa)':>12} {'P_tot(Pa)':>12} {'T_tot(K)':>9} "
          f"{'R*Vt':>8} {'rho':>8}")
    print("-"*100)

    for j in range(JM):
        vx  = float(np.mean(Vx[:,  j, K_mid]))
        vm  = float(np.mean(Vm[:,  j, K_mid]))
        vt  = float(np.mean(Vt[:,  j, K_mid]))
        ma  = float(np.mean(Mach_abs[:, j, K_mid]))
        ps  = float(np.mean(P_static[:, j, K_mid]))
        pt  = float(np.mean(P_total[:,  j, K_mid]))
        tt  = float(np.mean(T_total[:,  j, K_mid]))
        rvt = float(np.mean(RVt[:,      j, K_mid]))
        rho = float(np.mean(DEN[:,      j, K_mid]))
        print(f"{j+1:>5} {vx:>8.2f} {vm:>8.2f} {vt:>8.2f} {ma:>7.4f} "
              f"{ps:>12.1f} {pt:>12.1f} {tt:>9.2f} "
              f"{rvt:>8.4f} {rho:>8.5f}")

    # Mass flow at each J - simple integration
    print(f"\n=== MASS FLOW AT EACH J STATION ===")
    print(f"{'J':>5} {'mdot (kg/s)':>14}")
    print("-"*22)
    for j in range(JM):
        # Integrate rho*Vx over annular area
        # Use trapezoidal rule in K direction, sum over I (pitchwise)
        # Area element at each (i,k): dr * d(rtheta)/r * r = dr * dtheta * r
        mdot = 0.0
        for k in range(KM-1):
            dr   = float(np.mean(np.abs(R[:, j, k+1] - R[:, j, k])))
            drt  = float(np.mean(np.abs(RT[:, j, k])))  # pitch = 2*pi*r/N_blades
            rvx  = float(np.mean(DENVX[:, j, k]))
            mdot += rvx * dr * drt
        print(f"{j+1:>5} {mdot:>14.5f}")

if __name__ == '__main__':
    import sys
    fname = sys.argv[1] if len(sys.argv) > 1 else 'flow_out'
    read_flow_out(fname)
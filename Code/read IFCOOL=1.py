import numpy as np
import struct


# Readers and setup

gamma = 1.333
Rgas = 287.0
Cp = 1148.86

RPM = 8283.0
omega = 2*np.pi*RPM/60
W41 = 12.31

IM = 37
JM = 434
KM = 37

def read_record(f):
    marker = f.read(4)
    if len(marker) < 4:
        return None

    nbytes = struct.unpack('i', marker)[0]

    data = f.read(nbytes)

    f.read(4)

    return data

def read_flow_out(filename='Other files/flow_out'):

    with open(filename, 'rb') as f:

        NSTEP = struct.unpack('i', read_record(f))[0]

        print("NSTEP =", NSTEP)

        def arr():

            rec = read_record(f)

            a = np.frombuffer(rec, dtype=np.float32)

            return a.reshape((IM, JM, KM), order='F')

        RO      = arr()
        ROVX    = arr()
        ROVR    = arr()
        ROVT    = arr()
        ROE     = arr()
        ROSUB   = arr()
        VISCVAR = arr()
        SPARE   = arr()

    return RO, ROVX, ROVR, ROVT, ROE

def read_grid_out(filename='Other files/grid_out'):

    with open(filename, 'rb') as f:

        read_record(f)

        IM2, JM2, KM2 = struct.unpack('iii', read_record(f))

        print(IM2, JM2, KM2)

        CPFILE, GAFILE = struct.unpack('ff', read_record(f))

        INDLETE = np.frombuffer(read_record(f), dtype=np.int32)

        WRAD = np.frombuffer(read_record(f), dtype=np.float32)

        NBLADE = np.frombuffer(read_record(f), dtype=np.int32)

        X = np.zeros((JM, KM))

        R = np.zeros((JM, KM))

        RTHETA = np.zeros((IM, JM, KM))

        for j in range(JM):

            for k in range(KM):

                a = np.frombuffer(read_record(f), dtype=np.float32)

                X[j, k] = a[0]

                R[j, k] = a[1]

                RTHETA[:, j, k] = a[2:]

    return X, R, RTHETA, NBLADE


def compute_mdot(R, rho, Vm, RTHETA, NBLADE, blade_row_boundaries):

    mdot_J = np.zeros(JM)

    # Build a per-J nblade array
    nblade_J = np.ones(JM, dtype=float) * NBLADE[0]  # fallback
    for j_start, j_end, nb in blade_row_boundaries:
        nblade_J[j_start:j_end+1] = nb

    for j in range(JM):
        mdot = 0.0
        for i in range(IM-1):
            for k in range(KM-1):
                rho_avg = 0.25*(rho[i,j,k] + rho[i+1,j,k] +
                                rho[i,j,k+1] + rho[i+1,j,k+1])
                Vm_avg  = 0.25*(Vm[i,j,k]  + Vm[i+1,j,k]  +
                                Vm[i,j,k+1]  + Vm[i+1,j,k+1])

                dr       = R[j, k+1] - R[j, k]          # spanwise
                dRtheta  = RTHETA[i+1,j,k] - RTHETA[i,j,k]  # pitchwise

                dA = dr * dRtheta
                mdot += rho_avg * Vm_avg * dA

        mdot_J[j] = mdot * nblade_J[j]   # scale to full annulus

    return mdot_J


def get_stageResults(blade_row_boundaries, Vx, Vt, Wt,
                     T_static, T_total,
                     P_static, P_total,
                     R2D):

    alpha = np.degrees(np.arctan2(Vt, Vx))
    beta  = np.degrees(np.arctan2(Wt, Vx))

    exp = (gamma - 1.0)/gamma

    stations = {
        "Hub":  0,
        "Mean": KM//2,
        "Tip":  KM-1
    }
    stage_results = []

    # Stage quantities

    for stage in range(2):

        stator = blade_row_boundaries[2*stage]
        rotor  = blade_row_boundaries[2*stage+1]

        J1   = stator[0] + 5
        J2   = stator[1] - 5
        Jmid = rotor[0] + 5
        J3   = rotor[1] - 5
        Jm   = (rotor[0] + rotor[1])//2

        for name, k in stations.items():

            alpha1 = np.mean(alpha[:,J1,k])
            alpha2 = np.mean(alpha[:,J2,k])

            beta2  = np.mean(beta[:,Jmid,k])
            beta3  = np.mean(beta[:,J3,k])

            alpha3 = np.mean(alpha[:,J3,k])
            

            U = omega*R2D[Jm,k]

            phi = np.mean(Vx[:,Jm,k])/U

            dh0 = Cp*(
                np.mean(T_total[:,Jmid,k]) -
                np.mean(T_total[:,J3,k])
            )

            psi = dh0/U**2

            # NASA reaction
            Pt0    = np.mean(P_total[:,J1,k])
            Ps_mid = np.mean(P_static[:,Jmid,k])
            Ps_out = np.mean(P_static[:,J3,k])

            numerator   = 1.0 - (Ps_mid/Pt0)**exp
            denominator = 1.0 - (Ps_out/Pt0)**exp

            reaction = 1.0 - numerator/denominator

            stage_results.append([ stage+1, name, alpha1, alpha2, beta2, beta3, alpha3, phi, psi, reaction])

    # Whole turbine quantities

    k = KM//2

    J_in  = blade_row_boundaries[0][0] + 5
    J_out = blade_row_boundaries[-1][1] - 5

    Pt4  = np.mean(P_total[:,J_in,k])
    Ps42 = np.mean(P_static[:,J_out,k])
    Tt41 = np.mean(T_total[:,J_in,k])
    Tt42 = np.mean(T_total[:,J_out,k])

    PR_tt = Pt4/np.mean(P_total[:,J_out,k])
    PR_ts = Pt4/Ps42

    # # Overall efficiency
    # Tt42s = Tt41*(1.0/PR_tt)**((gamma-1)/gamma)
    # eta = (Tt41-Tt42)/(Tt41-Tt42s)

    # Blade speed ratio
    U = omega*np.mean(R2D[:,k])

    C0 = np.sqrt(
        2.0*Cp*Tt41*
        (1.0-(Ps42/Pt4)**exp)
    )

    U_C0 = U/C0

    # Overall loading
    psi_total = Cp*(Tt41-Tt42)/U**2

    T4  = np.mean(T_static[:, blade_row_boundaries[0][0]+15, k])
    T41 = np.mean(T_static[:, blade_row_boundaries[1][0]+5, k])
    T42 = np.mean(T_static[:, blade_row_boundaries[-1][1]-15, k])

    EnergyExtraction = Cp * (T4 - T42) / Tt41

    CorrectedSpeed = (RPM*2*np.pi/60) / np.sqrt(Tt41)

    Pt4 = np.mean(P_total[:, blade_row_boundaries[0][0]+5, k])
    FlowFunction = W41*np.sqrt(Tt41)/(Pt4/1000.0)


    turbine_results = {
        # "Efficiency": eta,
        "U/C0": U_C0,
        "Loading": psi_total,
        "PR_tt": PR_tt,
        "PR_ts": PR_ts,

        "EnergyExtraction": EnergyExtraction,
        "CorrectedSpeed": CorrectedSpeed,
        "FlowFunction": FlowFunction,
    }

    return stage_results, turbine_results


RO, ROVX, ROVR, ROVT, ROE = read_flow_out("Other files/flow_out_design1")
X2D, R2D, RTHETA, NBLADE = read_grid_out("Other files/grid_out_design1")
R = np.repeat(R2D[np.newaxis, :, :], IM, axis=0)

blade_row_boundaries = [
    (0,   110, 46),
    (111, 221, 76),
    (222, 310, 48),
    (311, 433, 70),
]


# Absolute quantities

rho = RO
rho_safe = np.where(np.abs(rho) < 1e-12, 1e-12, rho)
Vx = ROVX / rho_safe
Vr = ROVR / rho_safe
Vt = ROVT / rho_safe
Vm = np.sqrt(Vx**2 + Vr**2)
Vabs = np.sqrt(Vx**2 + Vr**2 + Vt**2)
E = ROE / rho_safe
T_static = (E - 0.5*Vabs**2)/(Cp - Rgas)
P_static = rho*Rgas*T_static
T_total = T_static + Vabs**2/(2*Cp)
P_total = P_static*(T_total/T_static)**(gamma/(gamma-1))


# Relative quantities

U = omega*R
Wt = Vt - U
Wrel = np.sqrt(Vx**2 + Vr**2 + Wt**2)
a = np.sqrt(np.maximum(gamma*Rgas*T_static, 0.0))
Mach_abs = Vabs/a
Mach_rel = Wrel/a
Swirl_angle = np.degrees(np.arctan2(Wt, Vm))
Pitch_angle = np.degrees(np.arctan2(Vr, Vm))


# MULTALL variables

Variable1 = 100*(Vm - np.mean(Vm[:,0,:]))/np.mean(Vm[:,0,:])
Variable2 = Vx
Variable3 = Vt
Variable4 = Vr
Variable5 = P_static
Variable6 = T_static
Variable7 = Mach_rel
Variable8 = T_total
Variable9 = Vm
Variable10 = Swirl_angle
Variable11 = Pitch_angle
Variable12 = rho
Variable13 = Mach_abs


# Plotting contours
import matplotlib.pyplot as plt

def plot_JK(variable, X2D, R2D, i_plane, title, cmap='jet'):
    """
    J-K plane at constant I
    """

    X = X2D
    R = R2D

    Z = variable[i_plane, :, :]

    plt.figure(figsize=(9,6))

    plt.contourf(X, R, Z, levels=40, cmap=cmap)

    plt.colorbar(label=title)

    plt.xlabel("Axial distance X (m)")
    plt.ylabel("Radius R (m)")
    plt.title(f"{title} (I = {i_plane+1})")

    plt.axis('equal')
    plt.tight_layout()
    plt.show()


def plot_IJ(variable, X2D, R2D, RTHETA, k_plane, title, cmap='jet'):
    """
    I-J plane at constant K
    """

    theta = np.zeros((IM, JM))

    for i in range(IM):
        theta[i,:] = RTHETA[i,:,k_plane] / R2D[:,k_plane]

    X = np.tile(X2D[:,k_plane], (IM,1))

    Y = theta

    Z = variable[:,:,k_plane]

    plt.figure(figsize=(12,5))

    plt.contourf(X, Y, Z, levels=40, cmap=cmap)

    plt.colorbar(label=title)

    plt.xlabel("Axial distance X (m)")
    plt.ylabel(r"$\theta$ (rad)")
    plt.title(f"{title} (K = {k_plane+1})")

    plt.tight_layout()
    plt.show()


variables = {
    "%Vm": Variable1,
    "Vx": Variable2,
    "Vtheta": Variable3,
    "Vr": Variable4,
    "Static Pressure": Variable5,
    "Relative Mach": Variable6,
    "Total Temperature": Variable7,
    "Meridional Velocity": Variable8,
    "Swirl Angle": Variable9,
    "Pitch Angle": Variable10,
    "Density": Variable11,
}


# J-K (annulus) contour
plot_JK(P_static, X2D, R2D, IM//2, "Static pressure")
# J-I (blade passage) contour
plot_IJ(P_static, X2D, R2D, RTHETA, KM//2, "Static pressure")

MassFlux = rho * Vm
# plot_JK(MassFlux, X2D, R2D, IM//2, "Mass Flux (kg/m²/s)")

# plot_IJ(MassFlux, X2D, R2D, RTHETA, KM//2, "Mass Flux (kg/m²/s)")

mdot_IJ = np.zeros((IM, JM, KM))

k_plane = KM // 2      # or hub/tip

for j in range(JM):

    mdot = 0.0

    for i in range(IM-1):

        rho_avg = 0.5 * (
            rho[i, j, k_plane] +
            rho[i+1, j, k_plane]
        )

        Vm_avg = 0.5 * (
            Vm[i, j, k_plane] +
            Vm[i+1, j, k_plane]
        )

        # Circumferential width
        dRtheta = (
            RTHETA[i+1, j, k_plane]
            - RTHETA[i, j, k_plane]
        )

        # Spanwise thickness (approximate using neighbouring K cells)
        if k_plane == 0:
            dr = R2D[j, 1] - R2D[j, 0]
        elif k_plane == KM-1:
            dr = R2D[j, KM-1] - R2D[j, KM-2]
        else:
            dr = 0.5 * (
                R2D[j, k_plane+1]
                - R2D[j, k_plane-1]
            )

        dA = dr * dRtheta

        mdot += rho_avg * Vm_avg * dA

        mdot_IJ[i+1, j] = mdot

plot_IJ(mdot_IJ, X2D, R2D, RTHETA, KM//2, "Cumulative Mass Flow (kg/s)")
# plt.figure(figsize=(8,6))
# plt.contourf(X2D, R2D, mdot_JK, 40)
# plt.colorbar(label="Integrated mass flow (kg/s)")
# plt.xlabel("X")
# plt.ylabel("R")
# plt.axis("equal")
# plt.show()


K_output = KM//2    # midspan
K_output = 0       # hub
K_output = KM-1    # tip

outfile = f"Output/Updated/IFCOOL = 1/multall_K{K_output+1}_COOL.txt"


with open(outfile, "w") as f:

    stage_results, turbine_results = get_stageResults(blade_row_boundaries, Vx, Vt, Wt, T_static, T_total, P_static, P_total, R2D)

    f.write("\nOVERALL TURBINE PERFORMANCE\n")
    f.write("="*60 + "\n")
    # f.write(f"Efficiency                {turbine_results['Efficiency']:.4f}\n")
    f.write(f"Blade-jet speed ratio     {turbine_results['U/C0']:.4f}\n")
    f.write(f"Pitchline loading         {turbine_results['Loading']:.4f}\n")
    f.write(f"Total-Total PR            {turbine_results['PR_tt']:.4f}\n")
    f.write(f"Total-Static PR           {turbine_results['PR_ts']:.4f}\n")
    f.write(f"Energy Extraction           {turbine_results['EnergyExtraction']:.4f}\n")
    f.write(f"Corrected Speed           {turbine_results['CorrectedSpeed']:.4f}\n")
    f.write(f"Flow function           {turbine_results['FlowFunction']:.4f}\n\n")

    f.write("STAGE PERFORMANCE\n")
    f.write("="*100 + "\n")
    f.write(
        "Stage  Location     Alpha1    Alpha2"
        "     Beta2     Beta3    Alpha3"
        "      Phi       Psi    Reaction\n"
    )

    for row in stage_results:

        f.write(
            f"{row[0]:5d}"
            f"{row[1]:>10}"
            f"{row[2]:11.2f}"   # alpha1
            f"{row[3]:11.2f}"   # alpha2
            f"{row[4]:11.2f}"   # alpha3
            f"{row[5]:11.2f}"   # beta2
            f"{row[6]:11.2f}"   # beta3
            f"{row[7]:10.4f}"   # phi
            f"{row[8]:10.4f}"   # psi
            f"{row[9]:12.4f}\n" # reaction
        )

    f.write("\n\n")
    
    f.write(f"Span station K = {K_output+1}\n\n")

    header = (
        "  J       %dVm        Vx(m/s)      Vtheta(m/s)    Vr(m/s)      "
        " Pstat(Pa)       Tstat(K)    Mrel      Mabs        T0(K)       "
        "Vm(m/s)     Swirl(deg)   Pitch(deg)  Density(kg/m3)     mdot(kg/s)\n"
    )

    f.write(header)
    f.write("-"*140 + "\n")

    Vm_ref = np.mean(Vm[:,0,:])

    for j in range(JM):

        pct_vm = np.mean(Variable1[:,j,K_output])

        vx = np.mean(Variable2[:,j,K_output])

        vt = np.mean(Variable3[:,j,K_output])

        vr = np.mean(Variable4[:,j,K_output])

        ps = np.mean(Variable5[:,j,K_output])

        ts = np.mean(Variable6[:,j,K_output])

        mrel = np.mean(Variable7[:,j,K_output])

        mabs = np.mean(Variable13[:,j,K_output])

        t0 = np.mean(Variable8[:,j,K_output])

        vm = np.mean(Variable9[:,j,K_output])

        swirl = np.mean(Variable10[:,j,K_output])

        pitch = np.mean(Variable11[:,j,K_output])

        rhoj = np.mean(Variable12[:,j,K_output])

        mdot = 0
        # print(j+1, end=" ")
        # mdot = compute_mdot(R2D, rho, Vm, RTHETA, NBLADE, blade_row_boundaries)[j]


        f.write(
            f"{j+1:4d}"
            f"{pct_vm:11.3f}"
            f"{vx:13.3f}"
            f"{vt:16.3f}"
            f"{vr:13.3f}"
            f"{ps:15.1f}"
            f"{ts:15.1f}"
            f"{mrel:10.4f}"
            f"{mabs:10.4f}"
            f"{t0:12.2f}"
            f"{vm:13.3f}"
            f"{swirl:13.3f}"
            f"{pitch:13.3f}"
            f"{rhoj:12.5f}"
            f"{mdot:19.3f}\n"
        )

print("Written to", outfile)
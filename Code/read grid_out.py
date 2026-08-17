import numpy as np
import struct

def read_record(f):
    marker = f.read(4)
    if len(marker) < 4:
        return None

    nbytes = struct.unpack('i', marker)[0]
    data = f.read(nbytes)
    f.read(4)

    return data


def read_grid_out(filename='grid_out'):

    with open(filename,'rb') as f:

        NSTEPS = struct.unpack('i', read_record(f))[0]

        IM,JM,KM = struct.unpack('iii', read_record(f))

        Cp,Ga = struct.unpack('ff', read_record(f))

        INDLETE = np.frombuffer(read_record(f),dtype=np.int32)
        WRAD     = np.frombuffer(read_record(f),dtype=np.float32)
        NBLADE   = np.frombuffer(read_record(f),dtype=np.int32)

        X = np.zeros((JM,KM))
        R = np.zeros((JM,KM))
        RTHETA = np.zeros((IM,JM,KM))

        for j in range(JM):
            for k in range(KM):

                arr = np.frombuffer(read_record(f),dtype=np.float32)

                X[j,k] = arr[0]
                R[j,k] = arr[1]

                RTHETA[:,j,k] = arr[2:]

    return X,R,RTHETA

X,R,RTHETA = read_grid_out("Other files/grid_out")

theta = np.zeros_like(RTHETA)
IM = 37
JM = 434
KM = 37

for j in range(JM):
    for k in range(KM):
        theta[:,j,k] = RTHETA[:,j,k]/R[j,k]

X3D = np.repeat(X[np.newaxis,:,:], IM, axis=0)
R3D = np.repeat(R[np.newaxis,:,:], IM, axis=0)

Y = R3D*np.cos(theta)
Z = R3D*np.sin(theta)



import matplotlib.pyplot as plt

# X and R are shape (JM,KM)

plt.figure(figsize=(12,5))

# constant K lines (hub → tip)
for k in range(KM):
    plt.plot(X[:,k], R[:,k], 'k-', linewidth=0.5)

# constant J lines (streamwise)
for j in range(JM):
    plt.plot(X[j,:], R[j,:], 'k-', linewidth=0.5)

plt.xlabel('X (m)')
plt.ylabel('R (m)')
plt.axis('equal')
plt.grid()
plt.show()

plt.figure(figsize=(12,4))

# constant k lines
k0 = 23

for i in range(IM):

    xx = []
    yy = []

    for j in range(JM):

        xx.append(X[j,k0])
        yy.append(RTHETA[i,j,k0])

    plt.plot(xx,yy,'k',linewidth=0.4)


# constant j lines
for j in range(JM):

    xx = []
    yy = []

    for i in range(IM):

        xx.append(X[j,k0])
        yy.append(RTHETA[i,j,k0])

    plt.plot(xx,yy,'k',linewidth=0.4)

plt.xlabel('X')
plt.ylabel('Rθ')
plt.axis('equal')
plt.show()
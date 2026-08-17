import numpy as np

INPUT_FILE = r"XRIN RIN source/XRIN,RIN values - final.txt"


# Read x,r pairs
data = np.loadtxt(INPUT_FILE, delimiter=",")

x = data[:, 0]
r = data[:, 1]


# Split regions
hub  = data[r < 33.0]
mean = data[(r >= 33.0) & (r <= 35.8)]
tip  = data[r > 35.8]


def write_region(filename, arr):

    # sort entire region by x
    arr = arr[np.argsort(arr[:, 0])]

    with open(filename, "w") as f:

        nblocks = len(arr) // 16

        for b in range(nblocks):

            block = arr[b*16:(b+1)*16]

            xvals = block[:, 0]/100
            rvals = block[:, 1]/100

            # x values
            f.write(" ".join(f"{v:11.6f}" for v in xvals[:8]) + "\n")
            f.write(" ".join(f"{v:11.6f}" for v in xvals[8:16]) + "\n")

            # r values
            f.write(" ".join(f"{v:11.6f}" for v in rvals[:8]) + "\n")
            f.write(" ".join(f"{v:11.6f}" for v in rvals[8:16]) + "\n")

            # LE / TE coordinates
            f.write(
                f"{xvals[2]:11.6f}"
                f"{xvals[13]:11.6f}"
                f"{rvals[2]:11.6f}"
                f"{rvals[13]:11.6f}"
                "  LEADING AND TRAILING EDGE COORDINATES\n"
            )


write_region("hub,mean,tip/Hub/hub.txt", hub)
write_region("hub,mean,tip/Mean/mean.txt", mean)
write_region("hub,mean,tip/Tip/tip.txt", tip)

print("hub.txt, mean.txt and tip.txt created.")
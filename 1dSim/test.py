import matplotlib.pyplot as plt
import glob

# Find all .dat files
files = glob.glob("result_*.dat")

for f in files:
    x, u = [], []
    itt, res = [], []
    residual = False
    with open(f) as file:
        for line in file:
            if line.startswith("RESIDUALS"):
                residual = True
                continue

            if residual:
                itti, resi = line.split()
                itt.append(int(itti))
                res.append(float(resi))

            else:
                xi, ui = line.split()
                x.append(float(xi))
                u.append(float(ui))

        residual = False
        fig, (ax1, ax2) = plt.subplots(1, 2)

        ax1.plot(x, u)
        ax1.set_xlabel("x")
        ax1.set_ylabel("u(x)")
        ax1.set_title(f"Final result of fitting {f}")

        ax2.plot(itt, res)
        ax2.set_xlabel("Itteration")
        ax2.set_ylabel("Residual")
        ax2.set_title("Residual progression for {f}")

        plt.show()

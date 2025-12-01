import subprocess
import shutil
import os


LEGACY_LOCATION = "../legacySolution/src/"
LEGACY_EXEC = "ionosphere.exe"
DATA_FOLDER = "../ionosphereSimulation/data/"

FILES = ["coords_out.txt",
         "jr_out.txt",
         "sig_out.txt",
         "dsig_out.txt",
         "pot_out.txt",
         "e_out.txt"]


def runProcesses():
    subprocess.run(["make", "-C", LEGACY_LOCATION],
                   stdout=subprocess.DEVNULL,
                   check=True)

    if not os.path.exists(LEGACY_LOCATION + LEGACY_EXEC):
        print("Compiled file not found")
        return 0

    subprocess.run(["./" + LEGACY_EXEC],
                   cwd=LEGACY_LOCATION,
                   check=True)

    print()
    for FILE in FILES:
        if os.path.exists(LEGACY_LOCATION + FILE):
            print("Found " + FILE)
            print(shutil.move(LEGACY_LOCATION + FILE, DATA_FOLDER + FILE))
        else:
            print("Failed to find " + FILE)


if __name__ == "__main__":
    runProcesses()

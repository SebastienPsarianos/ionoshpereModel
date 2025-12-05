import subprocess
import shutil
import os
import glob

LEGACY_LOCATION = "../legacySolution/src/"
LEGACY_EXEC = "ionosphere.exe"
DATA_FOLDER = "../ionosphereSimulation/data/"
NEW_BUILD_LOCATION = "../ionosphereSimulation/build/"
NEW_EXEC = "IonosphereSolver"

FILES = ["sim_input.txt",
         "sig_out.txt",
         "dsig_out.txt",
         "pot_out.txt",
         "e_out.txt"]


def cleanup():
    txt_files = glob.glob(os.path.join(DATA_FOLDER, '*.txt'))

    for file_path in txt_files:
        try:
            os.remove(file_path)
            print(f"Deleted: {file_path}")
        except OSError as e:
            print(f"Error deleting {file_path}: {e}")


def runProcesses():
    subprocess.run(["make"],
                   cwd=LEGACY_LOCATION,
                   stdout=subprocess.DEVNULL,
                   check=True)

    if not os.path.exists(LEGACY_LOCATION + LEGACY_EXEC):
        print("Compiled file not found")
        return 0

    print("Compilation Successful")

    subprocess.run(["./" + LEGACY_EXEC],
                   cwd=LEGACY_LOCATION,
                   stdout=subprocess.DEVNULL,
                   check=True)

    print("Ran legacy simulation")

    for FILE in FILES:
        if os.path.exists(LEGACY_LOCATION + FILE):
            print("Found " + FILE + " and moved to" +
                  shutil.move(LEGACY_LOCATION + FILE, DATA_FOLDER + FILE))
        else:
            print("Failed to find " + FILE)

    subprocess.run(["cmake", "../"],
                   cwd=NEW_BUILD_LOCATION,
                   stdout=subprocess.DEVNULL,
                   check=True)

    subprocess.run(["make"],
                   cwd=NEW_BUILD_LOCATION,
                   stdout=subprocess.DEVNULL,
                   check=True)

    print("New simulation build successful")

    if not os.path.exists(LEGACY_LOCATION + LEGACY_EXEC):
        print("Compiled file not found")
        return 0

    subprocess.run(["./" + NEW_EXEC, "../data/sim_input.txt"],
                   cwd=NEW_BUILD_LOCATION,
                   check=True)


if __name__ == "__main__":
    cleanup()
    runProcesses()

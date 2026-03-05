import subprocess
import shutil
import os
import glob
import argparse

LEGACY_DIR = "../legacySolution/src/"
LEGACY_EXEC = "ionosphere.exe"
NEW_DIR = "../ionosphereSimulation/"
NEW_DATA_DIR = "../ionosphereSimulation/data/"
DATA_DIR = "../data"

LEGACY_OUTPUT = ["sim_input.txt", "sig_out.txt",
                 "dsig_out.txt", "pot_out.txt", "e_out.txt"]


def run(cmd, **kwargs):
    subprocess.run(cmd, check=True, **kwargs)


def cleanup():
    for f in glob.glob(os.path.join(NEW_DATA_DIR, "*.txt")):
        os.remove(f)


def run_pipeline(output_flag="-t"):
    run(["make", "ionosphere"], cwd=LEGACY_DIR, stdout=subprocess.DEVNULL)
    print("Legacy build successful")

    run(["./" + LEGACY_EXEC], cwd=LEGACY_DIR, stdout=subprocess.DEVNULL)
    print("Legacy simulation complete")

    for f in LEGACY_OUTPUT:
        src = LEGACY_DIR + f
        if os.path.exists(src):
            shutil.move(src, NEW_DATA_DIR + f)

    run(["./build.sh"], cwd=NEW_DIR, stdout=subprocess.DEVNULL)
    print("New simulation build successful")

    run(["./run.sh", output_flag, "../data/sim_input.txt"], cwd=NEW_DIR)

    os.makedirs(DATA_DIR, exist_ok=True)
    for src in glob.glob(NEW_DATA_DIR + "*"):
        shutil.move(src, os.path.join(DATA_DIR, os.path.basename(src)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run ionosphere simulation pipeline.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-t", "--tecplot", dest="output_flag",
                       action="store_const",
                       const="-t",
                       help="Output in TecPlot format (default)")
    group.add_argument("-m", "--matplotlib", dest="output_flag",
                       action="store_const",
                       const="-m",
                       help="Output in matplotlib format for plot2.py")
    parser.set_defaults(output_flag="-t")
    args = parser.parse_args()

    cleanup()
    run_pipeline(args.output_flag)

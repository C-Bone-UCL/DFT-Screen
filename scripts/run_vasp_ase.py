#!/usr/bin/env python3
import os, sys, shutil, argparse, glob, subprocess
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun

def check_finished(f="OUTCAR"):
    if not os.path.isfile(f):
        return False
    with open(f) as h:
        return any("reached required accuracy" in ln for ln in h)

def get_cycle():
    outs = sorted(glob.glob("Cycle*ISIF3.OUTCAR"))
    return 0 if not outs else int(outs[-1].split('_')[0][5:]) + 1

def run_step(structure, incar_settings, tag):
    print(f"\n--- Preparing VASP step: '{tag}' ---")

    # Use Pymatgen's MPRelaxSet to generate a complete and validated
    # set of VASP input files (INCAR, KPOINTS, POSCAR, POTCAR).
    # It automatically uses the VASP_PP_PATH environment variable for potentials.
    # The 'force_gamma=False' argument honors your original setting.
    calc_set = MPRelaxSet(
        structure,
        user_incar_settings=incar_settings,
        force_gamma=False,
        user_potcar_functional="PBE_54"
    )
    calc_set.write_input('.')
    
    print(f"--- Input files for {tag} written successfully. ---")

    # print first few lines of INCAR for debugging
    with open("INCAR", "r") as incar_file:
        print("--- INCAR contents: ---")
        for line in incar_file.readlines()[:10]:
            print(line.strip())
        print("--- End of INCAR contents ---")
    # print first few lines of KPOINTS for debugging
    with open("KPOINTS", "r") as kpoints_file:
        print("--- KPOINTS contents: ---")
        for line in kpoints_file.readlines()[:10]:
            print(line.strip())
        print("--- End of KPOINTS contents ---")
    # print first few lines of POSCAR for debugging
    with open("POSCAR", "r") as poscar_file:
        print("--- POSCAR contents: ---")
        for line in poscar_file.readlines()[:10]:
            print(line.strip())
        print("--- End of POSCAR contents ---")


    # Manually execute VASP using the command from the environment
    vasp_command_str = os.environ["VASP_COMMAND"]
    command_list = vasp_command_str.split()
    print(f"--- Executing command: {' '.join(command_list)} ---")

    with open('vasp_out', 'w') as f_out:
        result = subprocess.run(command_list, stdout=f_out, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"--- VASP CRASHED for tag {tag} ---")
        print("--- stderr from VASP: ---")
        print(result.stderr)
        sys.exit(f"VASP execution failed for tag {tag}. Check outputs.")

    shutil.copy("OUTCAR",  f"{tag}.OUTCAR")
    shutil.copy("CONTCAR", f"{tag}.CONTCAR")
    shutil.copy("vasprun.xml", f"{tag}.vasprun.xml")

    # print first few lines of CONTCAR for debugging
    with open("CONTCAR", "r") as contcar_file:
        print("--- CONTCAR contents: ---")
        for line in contcar_file.readlines()[:10]:
            print(line.strip())
        print("--- End of CONTCAR contents ---")
    
    return "CONTCAR"

def workflow(cif, potcar_dir):
    # Read the initial structure using pymatgen
    structure = Structure.from_file(cif)
    structure.comment = f"Structure {os.path.splitext(cif)[0]}"

    # Manual POSCAR/POTCAR/CONTCAR creation is no longer needed.
    # Pymatgen handles all input file generation within each run_step call.

    ranks  = int(os.environ.get("NSLOTS", "1"))
    kpar   = max(1, int(round(ranks ** 0.5)))
    while ranks % kpar:
        kpar -= 1
    npar   = max(1, ranks // kpar)
    print(f"Using npar={npar} and kpar={kpar} for this run.")

    # Define INCAR settings as dictionaries. Pymatgen uses uppercase tags.
    common_settings = dict(
        PREC="Accurate",
        ENCUT=400,
        EDIFF=1e-6,
        EDIFFG=1e-5,
        ISMEAR=0,
        IBRION=2,
        ISPIN=1,
        LREAL=False,
        LWAVE=True,
        LCHARG=True,
        NELM=120,
        ISYM=2,
        SYMPREC=1e-8,
        NPAR=npar,
        KPAR=kpar,
    )

    isif2_settings  = {**common_settings, "ISIF": 2, "NSW": 60, "POTIM": 0.5}
    isif3_settings  = {**common_settings, "ISIF": 3, "NSW": 80, "POTIM": 0.75}
    isif3s_settings = {**common_settings, "ISIF": 3, "NSW": 8,  "POTIM": 0.75}

    cyc = get_cycle()
    while cyc < 20:
        # Run ISIF=2 relaxation (relax ions)
        run_step(structure, isif2_settings, f"Cycle{cyc}_ISIF2")
        structure = Structure.from_file("CONTCAR")

        # Run ISIF=3 relaxation (relax ions, cell shape, volume)
        run_step(structure, isif3_settings, f"Cycle{cyc}_ISIF3")
        structure = Structure.from_file("CONTCAR")

        if check_finished(f"Cycle{cyc}_ISIF3.OUTCAR"):
            # Final short relaxation to ensure convergence
            run_step(structure, isif3s_settings, f"Cycle{cyc}_ISIF3s")
            if check_finished(f"Cycle{cyc}_ISIF3s.OUTCAR"):
                print("--- Workflow converged successfully. ---")
                break
        cyc += 1
    
    # Analyze the final vasprun.xml from the last successful step
    vr = Vasprun("vasprun.xml", parse_eigen=True)
    gap = vr.get_band_structure().get_band_gap()["energy"]
    with open("results.txt", "w") as f:
        f.write(f"Energy_eV {vr.final_energy:.6f}\nBandGap_eV {gap:.4f}\n")
    os.chdir("..")

def main():
    pr = argparse.ArgumentParser()
    pr.add_argument("cif_dir")
    args = pr.parse_args()

    potdir = os.environ.get("VASP_PP_PATH")
    if not potdir or not os.environ.get("VASP_COMMAND"):
        sys.exit("export VASP_COMMAND and VASP_PP_PATH before running.")

    for cif in sorted(f for f in os.listdir(args.cif_dir) if f.lower().endswith(".cif")):
        case = os.path.splitext(cif)[0]
        if not os.path.isdir(case):
            os.mkdir(case)
        shutil.copy(os.path.join(args.cif_dir, cif), case)
        os.chdir(case)
        workflow(cif, potdir)
        os.chdir("..")

if __name__ == "__main__":
    main()
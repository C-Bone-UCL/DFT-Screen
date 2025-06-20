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

    calc_set = MPRelaxSet(
        structure,
        user_incar_settings=incar_settings,
        force_gamma=False,
        user_potcar_settings={"Ti": "Ti", "O": "O"}
    )
    
    calc_set.incar.write_file("INCAR")
    calc_set.kpoints.write_file("KPOINTS")
    structure.to(fmt="poscar", filename="POSCAR")
    
    potcar_dir = os.environ["PMG_VASP_PSP_DIR"]
    symbols_in_order = []
    for site in structure:
        symbol = site.specie.symbol
        if symbol not in symbols_in_order:
            symbols_in_order.append(symbol)
            
    with open("POTCAR", 'wb') as potcar_file:
        for symbol in symbols_in_order:
            potcar_path = os.path.join(potcar_dir, symbol, 'POTCAR')
            with open(potcar_path, 'rb') as individual_potcar:
                shutil.copyfileobj(individual_potcar, potcar_file)

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
    
    return "CONTCAR"

def workflow(cif, potcar_dir):
    structure = Structure.from_file(cif)
    structure.comment = f"Structure {os.path.splitext(cif)[0]}"

    ranks  = int(os.environ.get("NSLOTS", "1"))
    kpar   = max(1, int(round(ranks ** 0.5)))
    while ranks % kpar:
        kpar -= 1
    npar   = max(1, ranks // kpar)
    print(f"Using npar={npar} and kpar={kpar} for this run.")

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
        run_step(structure, isif2_settings, f"Cycle{cyc}_ISIF2")
        structure = Structure.from_file("CONTCAR")

        run_step(structure, isif3_settings, f"Cycle{cyc}_ISIF3")
        structure = Structure.from_file("CONTCAR")

        if check_finished(f"Cycle{cyc}_ISIF3.OUTCAR"):
            run_step(structure, isif3s_settings, f"Cycle{cyc}_ISIF3s")
            if check_finished(f"Cycle{cyc}_ISIF3s.OUTCAR"):
                print("--- Workflow converged successfully. ---")
                break
        cyc += 1
    
    vr = Vasprun("vasprun.xml", parse_eigen=True)
    gap = vr.get_band_structure().get_band_gap()["energy"]
    with open("results.txt", "w") as f:
        f.write(f"Energy_eV {vr.final_energy:.6f}\nBandGap_eV {gap:.4f}\n")
    # The os.chdir("..") line that was here has been removed.

def main():
    pr = argparse.ArgumentParser()
    pr.add_argument("cif_dir")
    args = pr.parse_args()

    potdir = os.environ.get("PMG_VASP_PSP_DIR")
    if not potdir or not os.environ.get("VASP_COMMAND"):
        sys.exit("export VASP_COMMAND and PMG_VASP_PSP_DIR before running.")

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
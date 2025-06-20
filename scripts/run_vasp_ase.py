#!/usr/bin/env python3
import os, sys, shutil, argparse, glob
from ase.io import read, write
from ase.calculators.vasp import Vasp
from pymatgen.io.vasp.outputs import Vasprun

def check_finished(f="OUTCAR"):
    if not os.path.isfile(f):
        return False
    with open(f) as h:
        return any("reached required accuracy" in ln for ln in h)

def get_cycle():
    outs = sorted(glob.glob("Cycle*_ISIF3.OUTCAR"))
    return 0 if not outs else int(outs[-1].split('_')[0][5:]) + 1

def run_step(calc, infile, tag):
    # --- Start Debug Block ---
    print(f"\n--- DEBUG: In run_step for tag '{tag}', about to read '{infile}' ---")
    if os.path.exists(infile):
        if os.path.getsize(infile) > 0:
            print(f"--- Contents of {infile}: ---")
            with open(infile, 'r') as f_debug:
                print(f_debug.read())
            print(f"--- End of {infile} contents ---")
        else:
            print(f"--- DEBUG WARNING: {infile} exists but is EMPTY. ---")
    else:
        print(f"--- DEBUG ERROR: {infile} does not exist! ---")
    # --- End Debug Block ---

    atoms = read(infile)
    atoms.calc = calc
    atoms.get_potential_energy()

    # --- Start Debug Block ---
    print(f"--- DEBUG: VASP run finished for tag '{tag}'. Checking outputs. ---")
    if os.path.exists("vasp_out"):
        print("--- DEBUG: Contents of vasp_out (VASP stdout) ---")
        with open("vasp_out", "r") as f_debug:
            # Print last 20 lines of VASP output
            lines = f_debug.readlines()
            for line in lines[-20:]:
                print(line.strip())
        print("--- END DEBUG: vasp_out ---")
    # --- End Debug Block ---

    shutil.copy("OUTCAR",  f"{tag}.OUTCAR")
    shutil.copy("CONTCAR", f"{tag}.CONTCAR")
    return "CONTCAR"

def workflow(cif, potcar_dir):
    atoms = read(cif)                             # already in your workflow
    atoms.info['comment'] = atoms.get_chemical_formula()   # <- clean first line
    # Write POSCAR. The vasp format writer sorts atoms by symbol by default.
    write("POSCAR", atoms, format="vasp", vasp5=True)

    # --- Start Manual POTCAR Generation ---
    # Get the exact order of elements from the written POSCAR
    symbols_in_order = []
    for atom in read("POSCAR"):
        symbol = atom.symbol
        if symbol not in symbols_in_order:
            symbols_in_order.append(symbol)

    # Concatenate POTCAR files in the correct order
    with open("POTCAR", 'wb') as potcar_file:
        for symbol in symbols_in_order:
            potcar_path = os.path.join(potcar_dir, symbol, 'POTCAR')
            with open(potcar_path, 'rb') as individual_potcar:
                shutil.copyfileobj(individual_potcar, potcar_file)
    # --- End Manual POTCAR Generation ---

    if not os.path.exists("CONTCAR"):
        shutil.copy("POSCAR", "CONTCAR")

    ranks  = int(os.environ.get("NSLOTS", "1"))
    kpar   = max(1, int(round(ranks ** 0.5)))
    while ranks % kpar:
        kpar -= 1
    npar   = max(1, ranks // kpar)

    print(f"Using npar={npar} and kpar={kpar} for this run.")

    common = dict(
    command=os.environ["VASP_COMMAND"],
    # We created the POTCAR manually, so we don't need pp, setups, or xc.
    # We specify the functional directly with gga='PE'.
    gga='PE',
    istart=0, icharg=2,
    prec="Accurate", lreal=False,
    encut=400, nelm=120, ediff=1e-6,
    kspacing=0.55, gamma=False,
    ispin=1, ediffg=1e-5, ibrion=2, isym=2, symprec=1e-8,
    ismear=0, lwave=True, lcharg=True,
    # --- FIX ---
    # Use the npar and kpar values calculated above
    npar=npar,
    kpar=kpar
    )

    isif2  = Vasp(**common, isif=2, nsw=60,  potim=0.5)
    isif3  = Vasp(**common, isif=3, nsw=80,  potim=0.75)
    isif3s = Vasp(**common, isif=3, nsw=8,   potim=0.75)

    cyc = get_cycle()
    while cyc < 20:
        run_step(isif2,  "CONTCAR", f"Cycle{cyc}_ISIF2")
        run_step(isif3,  "CONTCAR", f"Cycle{cyc}_ISIF3")
        if check_finished():
            run_step(isif3s, "CONTCAR", f"Cycle{cyc}_ISIF3s")
            if check_finished():
                break
        cyc += 1

    vr   = Vasprun("vasprun.xml", parse_eigen=True)
    gap  = vr.get_band_structure().get_band_gap()["energy"]
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
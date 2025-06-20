#!/usr/bin/env python3
import os, sys, shutil, argparse, glob, subprocess
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
    # Read the structure to be run
    atoms = read(infile)
    # Assign the calculator to the atoms object. This does NOT run the calculation.
    # It only sets up the parameters for writing the input files.
    atoms.calc = calc

    # --- Manually write VASP input files using ASE's writers ---
    # The calculator will write INCAR and KPOINTS. POSCAR and POTCAR are already present.
    atoms.calc.write_input(atoms)

    # --- Manually execute VASP ---
    # We use subprocess.run to have more control and capture output.
    vasp_command_str = os.environ["VASP_COMMAND"]
    # We split the command string into a list for subprocess
    # Note: This is a simple split; works for "mpirun -np 16 ..."
    # but might fail for more complex commands with quoted arguments.
    command_list = vasp_command_str.split()
    print(f"--- DEBUG: Executing command: {' '.join(command_list)} ---")

    # We redirect stdout to 'vasp_out' manually
    with open('vasp_out', 'w') as f_out:
        result = subprocess.run(command_list, stdout=f_out, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"--- VASP CRASHED for tag {tag} ---")
        print("--- stderr from VASP: ---")
        print(result.stderr)
        # We can decide to stop the workflow here if VASP crashes
        sys.exit(f"VASP execution failed for tag {tag}. Check outputs.")

    shutil.copy("OUTCAR",  f"{tag}.OUTCAR")
    shutil.copy("CONTCAR", f"{tag}.CONTCAR")
    return "CONTCAR"

def workflow(cif, potcar_dir):
    atoms = read(cif)                             # already in your workflow
    # This becomes the comment line in the POSCAR. Let's make it clean.
    atoms.info['comment'] = f"Structure {os.path.splitext(cif)[0]}"
    write("POSCAR", atoms, format="vasp", vasp5=True)

    # --- Start Manual POTCAR Generation ---
    symbols_in_order = []
    # Read the newly written POSCAR to get the exact atom order
    for atom in read("POSCAR"):
        symbol = atom.symbol
        if symbol not in symbols_in_order:
            symbols_in_order.append(symbol)

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
    # REMOVED command parameter as we now run VASP manually.
    # We tell ASE a POTCAR exists by not providing pp, setups, or xc.
    # Set GGA to PS for the PBEsol functional.
    gga='PS',
    istart=0, icharg=2,
    prec="Accurate", lreal=False,
    encut=400, nelm=120, ediff=1e-6,
    kspacing=0.55, gamma=False,
    ispin=1, ediffg=1e-5, ibrion=2, isym=2, symprec=1e-8,
    ismear=0, lwave=True, lcharg=True,
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
#!/usr/bin/env python3
import os, sys, shutil, argparse, glob
from ase.io import read, write
from ase.calculators.vasp import Vasp
from pymatgen.io.vasp.outputs import Vasprun

def check_finished(outfile="OUTCAR"):
    if not os.path.isfile(outfile):
        return False
    with open(outfile) as f:
        for line in f:
            if ("reached required accuracy"
                    in line and "energy minimisation" in line):
                return True
    return False

def get_cycle_index():
    outs = sorted(glob.glob("Cycle*_ISIF3.OUTCAR"))
    if not outs:
        return 0
    return int(outs[-1].split('_')[0][5:]) + 1

def run_one_step(calculator, in_file, tag):
    bulk = read(in_file)
    bulk.calc = calculator
    _ = bulk.get_potential_energy()
    shutil.copy("OUTCAR", f"{tag}.OUTCAR")
    shutil.copy("CONTCAR", f"{tag}.CONTCAR")
    return read("CONTCAR")

def workflow_for_cif(cif_file, vasp_exe, pp_path):
    os.makedirs("run", exist_ok=True)
    os.chdir("run")

    atoms = read(cif_file)
    write("POSCAR", atoms, format="vasp")
    if not os.path.exists("POSCAR_initial"):
        shutil.copy("POSCAR", "POSCAR_initial")
    if not os.path.exists("CONTCAR"):
        shutil.copy("POSCAR", "CONTCAR")

    os.environ["VASP_PP_PATH"] = pp_path
    command = f"mpirun -np 16 {vasp_exe} > vasp.log 2>&1"

    common = dict(
        command=command, istart=1, icharg=1, prec="Accurate", lreal=False,
        encut=400, nelm=120, ediff=1e-6, xc="PBE", kspacing=0.55,
        ispin=1, ediffg=1e-5, ibrion=2, isym=2, symprec=1e-8,
        ismear=0, lwave=True, lcharg=True, npar=8, kpar=2
    )

    isif2  = Vasp(**common, isif=2, nsw=150, potim=0.5)
    isif3  = Vasp(**common, isif=3, nsw=200, potim=0.75)
    isif3s = Vasp(**common, isif=3, nsw=8,   potim=0.75)

    cycle = get_cycle_index()
    while cycle < 30:  # hard stop
        atoms = run_one_step(isif2,  "CONTCAR", f"Cycle{cycle}_ISIF2")
        atoms = run_one_step(isif3,  "CONTCAR", f"Cycle{cycle}_ISIF3")

        if check_finished("OUTCAR"):
            atoms = run_one_step(isif3s, "CONTCAR", f"Cycle{cycle}_ISIF3s")
            if check_finished("OUTCAR"):
                break
        cycle += 1

    vr = Vasprun("vasprun.xml", parse_eigen=True)
    gap = vr.get_band_structure().get_band_gap()["energy"]
    with open("results.txt", "w") as f:
        f.write(f"Energy_eV {vr.final_energy:.6f}\nBandGap_eV {gap:.4f}\n")
    os.chdir("..")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("cif_folder")
    p.add_argument("--vasp_exe", required=True)
    p.add_argument("--potcar_dir", required=True)
    args = p.parse_args()

    cif_list = sorted([c for c in os.listdir(args.cif_folder) if c.lower().endswith(".cif")])
    for cif in cif_list:
        case = os.path.splitext(cif)[0]
        if not os.path.isdir(case):
            os.mkdir(case)
        shutil.copy(os.path.join(args.cif_folder, cif), case)
        os.chdir(case)
        workflow_for_cif(cif, args.vasp_exe, args.potcar_dir)
        os.chdir("..")

if __name__ == "__main__":
    main()

# Inspired from ASE for VASP in MDIG Group Guide

#!/usr/bin/env python3
import os, sys, glob, shutil, argparse
import ase
from ase.io import read, write
from ase.calculators.vasp import Vasp
from pymatgen.io.vasp.outputs import Vasprun

def get_max_ind(pattern):
    files = glob.glob(pattern)
    if not files:
        return 0
    return max(int(f.split('.')[1]) for f in files)

def check_finished(outfile):
    with open(outfile) as f:
        for line in f:
            if 'reached required accuracy - stopping structural energy minimisation' in line:
                return True
    return False

def run_workflow(cif_path, vasp_command, pp_path):
    atoms = ase.io.read(cif_path)

    os.makedirs("run", exist_ok=True)
    os.chdir("run")

    write("POSCAR", atoms, format="vasp")

    #### Set parameters here
    cutoff = 500.0 # This is the planewave cutoff energy in eV
    kp, np = 14, 6 # Number of k-points and processors
    max2, max4, max3, max3s = 2000, 10, 100, 10 # Max steps for each run
    sym, spin, ivdw = 2, 1, 0 # Symmetry, spin polarization, and vdW correction
    pim2, pim4, pim3 = 0.5, 0.05, 0.75 # Ionic relaxation parameters
    ioncut = 1e-6 # Ionic convergence criterion
    gga = 'pbesol' # Exchange-correlation functional
    symprec = 0.01 # Symmetry precision
    brion, ist = 1, 1 # Brion algorithm and initial structure
    #### Finished setting parameters

    #### Define the types of run that will be required
    isif2 = Vasp(system="System", istart=ist, iniwav=1, icharg=0, gamma=True, reciprocal=True,
                 prec="Accurate", lreal=False, algo="All", encut=cutoff,
                 nelm=200, ediff=1e-7, xc=gga, kspacing=0.242, ispin=spin,
                 ediffg=ioncut, nsw=max2, ibrion=brion, isif=2, isym=sym, ivdw=ivdw,
                 npar=np, kpar=kp, potim=pim2, symprec=symprec, ismear=0,
                 command=vasp_command, pp=pp_path, directory=".")

    isif3 = Vasp(system="System", istart=ist, iniwav=1, icharg=0, gamma=True, reciprocal=True,
                 prec="Accurate", lreal=False, algo="All", encut=cutoff,
                 nelm=200, ediff=1e-7, xc=gga, kspacing=0.242, ispin=spin,
                 ediffg=ioncut, nsw=max3, ibrion=brion, isif=3, isym=sym, ivdw=ivdw,
                 npar=np, kpar=kp, potim=pim3, symprec=symprec, ismear=0,
                 command=vasp_command, pp=pp_path, directory=".")

    isif3s = Vasp(system="System", istart=ist, iniwav=1, icharg=0, gamma=True, reciprocal=True,
                  prec="Accurate", lreal=False, algo="All", encut=cutoff,
                  nelm=200, ediff=1e-7, xc=gga, kspacing=0.242, ispin=spin,
                  ediffg=ioncut, nsw=max3s, ibrion=brion, isif=3, isym=sym, ivdw=ivdw,
                  npar=np, kpar=kp, potim=pim3, symprec=symprec, ismear=0,
                  command=vasp_command, pp=pp_path, directory=".")


    bulk = read("POSCAR")
    finished = False
    counter = get_max_ind('OUTCAR.*.i3')

    while not finished:
        bulk.calc = isif2
        bulk.get_potential_energy()
        shutil.copy("OUTCAR", f"OUTCAR.{counter}.i2")
        shutil.copy("CONTCAR", f"CONTCAR.{counter}.i2")
        bulk = read("CONTCAR")

        bulk.calc = isif3
        bulk.get_potential_energy()
        shutil.copy("OUTCAR", f"OUTCAR.{counter}.i3")
        shutil.copy("CONTCAR", f"CONTCAR.{counter}.i3")
        bulk = read("CONTCAR")

        bulk.calc = isif3s
        bulk.get_potential_energy()
        shutil.copy("OUTCAR", f"OUTCAR.{counter}.i3s")
        shutil.copy("CONTCAR", f"CONTCAR.{counter}.i3s")

        finished = check_finished("OUTCAR")
        bulk = read("CONTCAR")
        counter += 1

    vr = Vasprun("vasprun.xml", parse_eigen=True)
    gap = vr.get_band_structure().get_band_gap()["energy"]
    energy = vr.final_energy
    with open("results.txt", "w") as f:
        f.write(f"Energy_eV {energy:.6f}\nBandGap_eV {gap:.4f}\n")
    os.chdir("..")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("cif_dir")
    args = p.parse_args()
    vasp_cmd = os.environ["VASP_COMMAND"]
    pp_path = os.environ["VASP_PP_PATH"]
    for cif in sorted(f for f in os.listdir(args.cif_dir) if f.lower().endswith(".cif")):
        print(f"Processing {cif}...")

        name = os.path.splitext(cif)[0]
        if not os.path.exists(name):
            os.mkdir(name)
        shutil.copy(os.path.join(args.cif_dir, cif), name)
        os.chdir(name)
        run_workflow(cif, vasp_cmd, pp_path)
        os.chdir("..")

if __name__ == "__main__":
    main()


#!/usr/bin/env python3
import os, sys, argparse
from ase.io import read, write
from ase.calculators.vasp import Vasp
from pymatgen.io.vasp.outputs import Vasprun

def run_one(cif_path, vasp_command, pp_path):
    base = os.path.splitext(os.path.basename(cif_path))[0]
    workdir = f"{base}_calc"
    os.makedirs(workdir, exist_ok=True)

    atoms = read(cif_path)
    write(os.path.join(workdir, "POSCAR"), atoms, format="vasp")

    calc = Vasp(
        command = vasp_command,
        directory = workdir,
        xc = "PBE",
        encut = 400,
        kpts = [2,2,2],
        istart = 0,
        icharg = 2,
        npar = 4,
        ibrion = -1,
        nsw = 0,
        ediff = 1e-5,
        lreal = "Auto",
        pp = pp_path
    )
    atoms.calc = calc

    print(f"Running VASP for {base} …")
    energy = atoms.get_potential_energy()

    vr = Vasprun(os.path.join(workdir, "vasprun.xml"), parse_dos=False, parse_eigen=True)
    try:
        gap = vr.get_band_structure().get_band_gap()["energy"]
    except Exception:
        gap = None

    with open(os.path.join(workdir, "results.txt"), "w") as f:
        f.write(f"Total energy (eV): {energy:.6f}\n")
        f.write(f"Indirect band gap (eV): {gap if gap is not None else 'N/A'}\n")

    write(os.path.join(workdir, "CONTCAR"), atoms, format="vasp")
    print(f"Finished {base}\n")

def main():
    p = argparse.ArgumentParser(description="VASP ground-state for all CIFs in a folder")
    p.add_argument("cif_dir", metavar="CIF_DIR")
    args = p.parse_args()

    vasp_cmd = os.environ.get("VASP_COMMAND")
    pp_path  = os.environ.get("VASP_PP_PATH")
    if not (vasp_cmd and pp_path):
        sys.exit("VASP_COMMAND and/or VASP_PP_PATH not set")

    if not os.path.isdir(args.cif_dir):
        sys.exit(f"{args.cif_dir} is not a directory")

    cif_files = [os.path.join(args.cif_dir, f)
                 for f in os.listdir(args.cif_dir) if f.lower().endswith(".cif")]

    if not cif_files:
        sys.exit("No CIF files found")

    for cif in sorted(cif_files):
        run_one(cif, vasp_cmd, pp_path)

if __name__ == "__main__":
    main()

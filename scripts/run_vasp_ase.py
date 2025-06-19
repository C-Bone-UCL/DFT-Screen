#!/usr/bin/env python3
import os
import sys
import argparse

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
        command=vasp_command,
        directory=workdir,
        xc="PBE",
        encut=400,
        kpts=[2,2,2],
        istart=0,
        icharg=2,
        npar=4,
        ibrion=-1,
        nsw=0,
        ediff=1e-5,
        lreal="Auto",
        pp=pp_path
    )
    atoms.set_calculator(calc)

    print(f"Running VASP for {base} …")
    energy = atoms.get_potential_energy()

    vr = Vasprun(os.path.join(workdir, "vasprun.xml"), parse_dos=False, parse_eigen=False)
    bg_data = vr.get_band_structure(line_mode=True).get_band_gap()
    band_gap = bg_data["energy"]

    out = os.path.join(workdir, "results.txt")
    with open(out, "w") as f:
        f.write(f"Total energy (eV): {energy:.6f}\n")
        f.write(f"Indirect band gap (eV): {band_gap:.4f}\n")

    # alse write out the structure
    write(os.path.join(workdir, "CONTCAR"), atoms, format="vasp")

    print(f"Results for {base} written to {out}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Run VASP ground-state & band-gap for all CIFs in a folder"
    )
    parser.add_argument(
        "cif_dir", metavar="CIF_DIR",
        help="Directory containing CIF files"
    )
    args = parser.parse_args()

    vasp_cmd = os.environ.get("VASP_COMMAND")
    pp_path  = os.environ.get("VASP_PP_PATH")
    if not vasp_cmd or not pp_path:
        print("Error: set VASP_COMMAND and VASP_PP_PATH in your env", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.cif_dir):
        print(f"Error: not a directory: {args.cif_dir}", file=sys.stderr)
        sys.exit(1)

    cif_files = [os.path.join(args.cif_dir, f)
                 for f in os.listdir(args.cif_dir)
                 if f.lower().endswith(".cif")]

    if not cif_files:
        print(f"No CIF files found in {args.cif_dir}", file=sys.stderr)
        sys.exit(0)

    for cif in sorted(cif_files):
        run_one(cif, vasp_cmd, pp_path)

if __name__ == "__main__":
    main()

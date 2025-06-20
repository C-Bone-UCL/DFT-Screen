#!/usr/bin/env python3
import os, sys, argparse, glob
from ase.io import read, write
from ase.calculators.vasp import Vasp
from pymatgen.io.vasp.outputs import Vasprun

def run_relaxation(cif_file):
    """
    Performs a full VASP structural relaxation for a single CIF file.
    """
    # Generate a clean name for the calculation directory
    base_name = os.path.splitext(os.path.basename(cif_file))[0]
    calc_dir = f"{base_name}_calc"

    print(f"\n--- Starting relaxation for: {base_name} in directory: {calc_dir} ---")

    # 1. Read the initial structure from the CIF file
    atoms = read(cif_file)

    # 2. Set up the VASP calculator for a full relaxation
    ranks = int(os.environ.get("NSLOTS", "1"))
    kpar = max(1, int(round(ranks**0.5)))
    while ranks % kpar:
        kpar -= 1
    npar = max(1, ranks // kpar)

    calc = Vasp(
        # --- File and Execution Control ---
        directory=calc_dir,
        command=os.environ.get("ASE_VASP_COMMAND"),
        txt=os.path.join(calc_dir, "vasp.out"),

        # --- Functional and Potential Selection (The "Legacy" Way) ---
        # We explicitly provide the path to the potential directory.
        # This overrides any automatic searching ASE would do.
        pp=os.environ.get("VASP_PP_PATH"),
        # Since we use 'pp', we should set the functional with 'gga' not 'xc'
        # to avoid conflicts. 'PS' is the VASP tag for PBEsol.
        gga='PS',

        # --- Accuracy ---
        prec='Accurate',
        encut=400,
        ediff=1e-6,

        # --- Relaxation Parameters ---
        ibrion=2,
        isif=3,
        nsw=80,
        ediffg=-0.01,

        # --- Electronic Step and k-points ---
        istart=0, icharg=2,
        nelm=120,
        kspacing=0.55,
        gamma=False,

        # --- Parallelization ---
        npar=npar,
        kpar=kpar,

        # --- Other Settings ---
        ispin=1,
        lreal=False,
        ismear=0,
        lwave=True,
        lcharg=True,
    )

    # 3. Attach the calculator to the atoms object
    atoms.calc = calc

    # 4. Run the calculation
    final_energy = atoms.get_potential_energy()
    print(f"--- Relaxation finished. Final Energy: {final_energy:.6f} eV ---")

    # 5. Post-process the results
    final_structure_file = os.path.join(calc_dir, f"{base_name}_relaxed.vasp")
    write(final_structure_file, atoms, format='vasp')
    print(f"Final relaxed structure saved to: {final_structure_file}")

    vasprun_path = os.path.join(calc_dir, "vasprun.xml")
    if os.path.exists(vasprun_path):
        vr = Vasprun(vasprun_path, parse_eigen=True)
        gap = vr.get_band_structure().get_band_gap()["energy"]
        print(f"Band Gap: {gap:.4f} eV")

        results_file = os.path.join(calc_dir, "results.txt")
        with open(results_file, "w") as f:
            f.write(f"Energy_eV {final_energy:.6f}\n")
            f.write(f"BandGap_eV {gap:.4f}\n")
        print(f"Summary saved to: {results_file}")
    else:
        print("--- WARNING: vasprun.xml not found. Cannot determine band gap. ---")


def main():
    parser = argparse.ArgumentParser(
        description="Run a VASP relaxation for all CIF files in a directory."
    )
    parser.add_argument("cif_dir", help="Directory containing input CIF files.")
    args = parser.parse_args()

    if not os.environ.get("VASP_PP_PATH") or not os.environ.get("ASE_VASP_COMMAND"):
        sys.exit("Error: Set VASP_PP_PATH and ASE_VASP_COMMAND environment variables.")

    cif_files = sorted(glob.glob(os.path.join(args.cif_dir, '*.cif')))
    if not cif_files:
        sys.exit(f"No .cif files found in '{args.cif_dir}'")

    for cif_file in cif_files:
        try:
            run_relaxation(cif_file)
        except Exception as e:
            print(f"\n--- An error occurred while processing {cif_file} ---")
            print(f"Error: {e}")
            print("--- Skipping to next file. ---")

if __name__ == "__main__":
    main()
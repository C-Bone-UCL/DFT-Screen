#!/usr/bin/env python3
import os, sys, shutil, argparse, glob, subprocess, time
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun

def check_finished(f="OUTCAR"):
    """
    Check if the VASP run has reached the required accuracy.
    """

    if not os.path.isfile(f):
        return False
    with open(f) as h:
        return any("reached required accuracy" in ln for ln in h)

def get_cycle():
    """ 
    Inspects the filesystem for previously completed cycle outputs
    By finding the highest-numbered cycle that has already run, resumes 
    from the next logical step, rather than starting over from cycle 0
    """

    outs = sorted(glob.glob("Cycle*ISIF3.OUTCAR"))
    return 0 if not outs else int(outs[-1].split('_')[0][5:]) + 1

def run_step(structure, incar_settings, tag):
    """
    Run a single VASP step with the given structure and settings.
    
    1. Prepares the VASP input files (INCAR, KPOINTS, POSCAR, POTCAR).
    Using MPRelaxSet for standard settings.
    2. Executes the VASP command defined in the environment variable VASP_COMMAND.
    3. Checks for successful completion and copies output files.
    4. Returns the name of the CONTCAR file for the next step.
    """
    
    print(f"\n####Preparing VASP step: '{tag}'####\n")

    calc_set = MPRelaxSet(
        structure,
        user_incar_settings=incar_settings,
        force_gamma=False,
        user_potcar_settings={"Ti": "Ti", "O": "O"} # Adjust as needed for your elements
    )
    # Default is the PBE fucntional
    
    # Make INCAR, KPOINTS, and POSCAR files normally
    calc_set.incar.write_file("INCAR")
    calc_set.kpoints.write_file("KPOINTS")
    structure.to(fmt="poscar", filename="POSCAR")
    
    potcar_dir = os.environ["PMG_VASP_PSP_DIR"]

    # This part collects all unique species in the structure
    # and writes a combined POTCAR file
    # Done manually because pymatgen's default POTCAR generation
    # does not support our potentials directory structure
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

    print(f"####Executing command: {' '.join(command_list)}####")

    with open('vasp_out', 'w') as f_out:
        result = subprocess.run(command_list, stdout=f_out, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"####VASP CRASHED for tag {tag}####")
        print("####stderr from VASP:####")
        print(result.stderr)
        sys.exit(f"VASP execution failed for tag {tag}. Check outputs.")

    # to record relaxation trajectory
    shutil.copy("OUTCAR",  f"{tag}.OUTCAR")
    shutil.copy("CONTCAR", f"{tag}.CONTCAR")
    shutil.copy("vasprun.xml", f"{tag}.vasprun.xml")
    
    return "CONTCAR"

def workflow(cif, potcar_dir):
    """
    Main workflow function to run VASP calculations on a given CIF file.

    1. Reads the CIF file to create a pymatgen Structure object.
    2. Sets up the VASP calculation parameters based on the environment.
    3. Runs the VASP calculations in a loop until convergence is achieved.
        a. Uses ISIF=2 for initial relaxation of ions while keeping cell fixed.
        b. Uses ISIF=3 for further relaxation of both ions and cell.
        c. Uses ISIF=3(s) with a smaller number of steps for final convergence.
    4. Outputs the final energy and band gap to a results file.
    """

    structure = Structure.from_file(cif)
    structure.comment = f"Structure {os.path.splitext(cif)[0]}"

    # This bit is for parallelization settings
    # It determines the number of processors to use based on the NSLOTS environment variable.
    ranks  = int(os.environ.get("NSLOTS", "1"))
    kpar   = max(1, int(round(ranks ** 0.5)))
    while ranks % kpar:
        kpar -= 1
    npar   = max(1, ranks // kpar)
    print(f"Using npar={npar} and kpar={kpar} for this run.")

    # Common settings for all VASP runs
    common_settings = dict(
        PREC="Accurate",    # Internal precision related params
        ENCUT=700,          # Energy cutoff for plane waves (size of basis set, scaling is expensive)
        EDIFF=1e-6,         # Energy convergence criterion
        EDIFFG=1e-5,        # Force convergence criterion
        ISMEAR=0,           # smearing applied to electronic states (0 = Gaussian smearing)
        SIGMA=0.1,          # Width of smearing function
        IBRION=2,           # Ionic relaxation algorithm (2 = conjugate gradient)
        ISPIN=1,            # Spin polarization (1 = non-spin-polarized, 2 = if magnetic)
        LREAL=False,        # Whether to use real-space projection for faster calculations
        LWAVE=True,         # Write WAVECAR file (useful for restarting calculations)
        LCHARG=True,        # Write CHGCAR file (useful for charge density analysis)
        NELM=120,           # Max steps for electronic convergence
        ISYM=2,             # Symmetry handling (2 = automatic symmetry detection)
        SYMPREC=1e-5,       # Precision for symmetry detection
        NPAR=npar,          # Number of processors for parallelization
        KPAR=kpar,          # K-point parallelization
    )

    # ISIF : Degrees of freedom for relaxation
    # ISIF=2: Relax ions, keep cell fixed
    # ISIF=3: Relax both ions and cell
    # ISIF=3s: Relax both ions and cell with fewer steps for final convergence

    # NSW : Number of ionic steps
    # POTIM : Ionic step size (time step for ionic dynamics)

    isif2_settings  = {**common_settings, "ISIF": 2, "NSW": 60, "POTIM": 0.5}
    isif3_settings  = {**common_settings, "ISIF": 3, "NSW": 80, "POTIM": 0.75}
    isif3s_settings = {**common_settings, "ISIF": 3, "NSW": 8,  "POTIM": 0.75}

    cyc = get_cycle()
    while cyc < 20:
        run_step(structure, isif2_settings, f"Cycle{cyc}_ISIF2")
        structure = Structure.from_file("CONTCAR")

        # print intermediate results
        vr = Vasprun("vasprun.xml", parse_eigen=True)
        print(f"ISIF2 - Energy: {vr.final_energy:.6f} eV")

        run_step(structure, isif3_settings, f"Cycle{cyc}_ISIF3")
        structure = Structure.from_file("CONTCAR")

        # Print intermediate results
        vr = Vasprun("vasprun.xml", parse_eigen=True)
        print(f"ISIF3 - Energy: {vr.final_energy:.6f} eV")

        if check_finished(f"Cycle{cyc}_ISIF3.OUTCAR"):
            run_step(structure, isif3s_settings, f"Cycle{cyc}_ISIF3s")
            if check_finished(f"Cycle{cyc}_ISIF3s.OUTCAR"):
                print("####Workflow converged successfully.####")
                # Print final results
                vr = Vasprun("vasprun.xml", parse_eigen=True)
                print(f"ISIF3s Energy: {vr.final_energy:.6f} eV")
                print(f"ISIF3s Band Gap: {vr.get_band_structure().get_band_gap()['energy']:.4f} eV")
                break
        cyc += 1
    
    vr = Vasprun("vasprun.xml", parse_eigen=True)
    gap = vr.get_band_structure().get_band_gap()["energy"]
    return vr.final_energy, gap

def main():
    pr = argparse.ArgumentParser()
    pr.add_argument("cif_dir")
    args = pr.parse_args()

    potdir = os.environ.get("PMG_VASP_PSP_DIR")
    if not potdir or not os.environ.get("VASP_COMMAND"):
        sys.exit("export VASP_COMMAND and PMG_VASP_PSP_DIR before running.")

    output_dir = "vasp_outputs"
    os.makedirs(output_dir, exist_ok=True)

    for cif in sorted(f for f in os.listdir(args.cif_dir) if f.lower().endswith(".cif")):
        case = os.path.splitext(cif)[0]
        case_path = os.path.join(output_dir, case)
        
        if not os.path.isdir(case_path):
            os.makedirs(case_path)
            
        shutil.copy(os.path.join(args.cif_dir, cif), case_path)
        
        original_dir = os.getcwd()
        os.chdir(case_path)
        
        start_time = time.time()
        final_energy, band_gap = workflow(cif, potdir)
        end_time = time.time()
        duration = end_time - start_time
        
        with open("results.txt", "w") as f:
            f.write(f"Energy_eV {final_energy:.6f}\n")
            f.write(f"BandGap_eV {band_gap:.4f}\n")
            f.write(f"Runtime_s {duration:.2f}\n")
            
        os.chdir(original_dir)

if __name__ == "__main__":
    main()
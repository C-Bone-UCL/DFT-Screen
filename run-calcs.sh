#!/bin/bash -l
source /etc/profile.d/modules.sh

#$ -N optimize_TiO2
#$ -l h_rt=12:00:00
#$ -pe mpi 36
#$ -l mem=1G

# Set the working directory
#$ -wd $HOME/Scratch/DFT-Screen

module purge
module load beta-modules
module load gcc-libs/10.2.0
module load compilers/intel/2022.2
module load mpi/intel/2021.6.0/intel

eval "$(conda shell.bash hook)"

export NSLOTS=${NSLOTS:=1}
export VASP_EXE=$HOME/Scratch/vasp/vasp-6.4.2/bin/vasp_std
export VASP_COMMAND="mpirun -np $NSLOTS $VASP_EXE"
export VASP_PP_PATH=$HOME/Scratch/vasp/vasp-6.4.2/potentials/potentials_PBE.54
export PMG_VASP_PSP_DIR=$VASP_PP_PATH

conda activate DFT-Screen-venv

python scripts/run_vasp_ase.py \
      files/TiO2_run
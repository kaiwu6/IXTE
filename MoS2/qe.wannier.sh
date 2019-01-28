#PBS -N wannier
#PBS -q debug
#PBS -l mppwidth=1
#PBS -l walltime=0:15:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load espresso/5.1.2
module load wannier90/2.0

aprun -n 1 wannier90.x MoS2

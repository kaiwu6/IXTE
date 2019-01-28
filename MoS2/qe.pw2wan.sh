#PBS -N pw2wan
#PBS -q debug
#PBS -l mppwidth=24
#PBS -l walltime=0:10:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load espresso/5.1.2
module load wannier90/2.0

aprun -n 24 pw2wannier90.x < MoS2.pw2wan > MoS2.pw2wan.out

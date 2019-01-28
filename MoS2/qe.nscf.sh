#PBS -N nscf
#PBS -q debug
#PBS -l mppwidth=24
#PBS -l walltime=0:30:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load espresso/5.1.2

aprun -n 24 pw.x < MoS2.nscf > MoS2.nscf.out

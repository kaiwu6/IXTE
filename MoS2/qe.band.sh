#PBS -N bands
#PBS -q debug
#PBS -l mppwidth=24
#PBS -l walltime=0:30:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load espresso/5.1.2

aprun -n 24 bands.x < bands.in > bands.out

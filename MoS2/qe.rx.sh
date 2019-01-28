#PBS -N relax
#PBS -q debug
#PBS -l mppwidth=24
#PBS -l walltime=0:30:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load espresso/5.1.2

#aprun -n 24 pw.x < MoS2_a600.rx > MoS2_a600.rx.out
#aprun -n 24 pw.x < MoS2_a601.rx > MoS2_a601.rx.out
#aprun -n 24 pw.x < MoS2_a602.rx > MoS2_a602.rx.out
#aprun -n 24 pw.x < MoS2_a603.rx > MoS2_a603.rx.out
#aprun -n 24 pw.x < MoS2_a603.rx > MoS2_a603.rx.out
aprun -n 24 pw.x < MoS2_a6035.rx > MoS2_a6035.rx.out
aprun -n 24 pw.x < MoS2_a605.rx > MoS2_a605.rx.out

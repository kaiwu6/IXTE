set style data dots
set nokey
set xrange [0: 3.18222]
set yrange [ -4.78199 :  2.98321]
set arrow from  0.56795,  -4.78199 to  0.56795,   2.98321 nohead
set arrow from  1.13590,  -4.78199 to  1.13590,   2.98321 nohead
set arrow from  1.78526,  -4.78199 to  1.78526,   2.98321 nohead
set arrow from  2.41484,  -4.78199 to  2.41484,   2.98321 nohead
set arrow from  3.08377,  -4.78199 to  3.08377,   2.98321 nohead
set arrow from  3.13299,  -4.78199 to  3.13299,   2.98321 nohead
set xtics (" G "  0.00000," S "  0.56795," M "  1.13590," K "  1.78526," L "  2.41484," G "  3.08377," D "  3.13299," A "  3.18222)
 plot "MoS2_band.dat"

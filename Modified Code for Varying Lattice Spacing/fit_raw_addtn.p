set xlabel "Lattice Spacing (s)"
set ylabel "Autocorrelation Time"

set border linewidth 1.5
set key right top
set terminal pngcairo
set output 'Fit_raw_addtn.dat.png'
plot 'raw_addtn.dat' t 'Data Points'

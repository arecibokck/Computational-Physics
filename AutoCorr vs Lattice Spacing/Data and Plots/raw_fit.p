set xlabel "Lattice Spacing (s)"
set ylabel "Autocorrelation Time"
a = -0.25
b = 0.4
func(x)=b*(x**a)
fit [:] func(x) 'raw.dat' via a,b

set border linewidth 1.5
set key right top
set terminal pngcairo
set output 'raw_fit.png'
plot 'raw.dat' with yerrorbars t 'Data Points', func(x) t 'Fit' 

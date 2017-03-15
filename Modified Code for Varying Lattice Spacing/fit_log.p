set xlabel "Log(Lattice Spacing)"
set ylabel "Log(Autocorrelation Time)"
a = -0.25
b = 0.5
func(x)=a*x+b
fit [:] func(x) 'log.dat' via a,b

set border linewidth 1.5
set terminal pngcairo
set output 'Fit_log.png'
plot 'log.dat' with yerrorbars t 'Data Points', func(x) t 'Fit' 

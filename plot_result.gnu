set term pngcairo

FINPUT = ARG1 # simulation data
FPLOT_PREFIX = 'res/plots/sim' # filename prefix for plots

# -------------------
# PLOTTING
# -------------------

# energy per lattice site
set output sprintf('%s_%s', FPLOT_PREFIX, 'e.png')

set title "Energy per lattice site"
set mxtics; set mytics
set nokey
set ylabel '<e>'; set xlabel 'T^*'
plot FINPUT u 1:2:3 w errorbars 

# susceptibility
set output sprintf('%s_%s', FPLOT_PREFIX, 'susc.png')

set title "Susceptibility"
set mxtics; set mytics
set nokey
set ylabel 'χ'; set xlabel 'T^*'
plot FINPUT u 1:8

# specific heat
set output sprintf('%s_%s', FPLOT_PREFIX, 'c.png')

set title "Specific heat"
set mxtics; set mytics
set nokey
set ylabel 'c^*'; set xlabel 'T^*'
plot FINPUT u 1:7 t 'c^*'
plot FINPUT u 1:9 t 'dE/dt'

# magnetization
set output sprintf('%s_%s', FPLOT_PREFIX, 'm.png')

set title "Specific heat"
set mxtics; set mytics
set nokey
set ylabel 'c^*'; set xlabel 'T^*'
plot FINPUT u 1:5:6 t '|m|' w yerrorbars 
plot FINPUT u 1:4 t 'sqrt(<m^2>)'
#!/usr/bin/env gnuplot



# Ising model



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced color size 9, 7 font "Helvetica, 36"
set output "is2_T2.eps"

reset

set style line 1 lt 1 lc rgb "#002080" lw 1.0 pt  7 ps 1.7
set style line 2 lt 2 lc rgb "#a00000" lw 1.5 pt  4 ps 1.4
set style line 3 lt 3 lc rgb "#a00040" lw 1.5 pt 12 ps 2.0
set style line 4 lt 4 lc rgb "#00a020" lw 1.5 pt  8 ps 2.0
set style line 5 lt 5 lc rgb "#a0a000" lw 1.5 pt 10 ps 2.0
set style line 7 lt 1 lc rgb "#808080"

set multiplot

set size 1, 1
set origin 0, 0

set xtics 1 offset 0, 0.2
set mxtics 2
set xlabel "Simulation time, {/Times-Italic t} ({/Symbol \264} 10^{/*0.7 3})" offset 0, 0.5
set xrange [0:10]

set ytics 1
set mytics 2
set ylabel 'Estimated entropy, ~{/Times-Italic S}{0.5\^}' offset 1.0, 0
set yrange [0:4.0]

set key left bottom Left reverse spacing 1.5 at 1.0, 1.7 width -5 maxrows 3

fn0 = "../../data/is2/is2_m0_T2.log"
fn1 = "../../data/is2/is2_m1_T2.log"
ref = 0
#print ref

plot [:][:] \
    fn0 u ($1/1e3): 2: 3 every 1 w lp ls 1 pt  6 ps 2   t "Metropolis, Uncorrected", \
    fn0 u ($1/1e3): 5: 6 every 1 w lp ls 2 pt  4 ps 1.4 t "Metropolis, Linear", \
    fn0 u ($1/1e3): 9:10 every 1 w lp ls 4 pt  8 ps 2   t "Metropolis, Exponential", \
    fn1 u ($1/1e3): 2: 3 every 1 w lp ls 1 pt  7 ps 2   t "Wolff, Uncorrected", \
    fn1 u ($1/1e3): 5: 6 every 1 w lp ls 2 pt  5 ps 1.4 t "Wolff, Linear", \
    fn1 u ($1/1e3): 9:10 every 1 w lp ls 4 pt  9 ps 2   t "Wolff, Exponential", \
    ref w l ls 7 t "Reference"


unset multiplot
unset output
set terminal pop
reset

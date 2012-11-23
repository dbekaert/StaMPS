# Note: different template for gnuplot 4.2, still i need to check
# with lower versions
#
# CHANGE.LOG: 
# * terminal png size defined
# * size of multiplot and plot changed 2->1,  1->0.5
# * zlabel manually set, since zlabel rotate doesnt work


set terminal png size 1024,768
set output "DUMPFILE" 

set size 1,1
set origin 0,0
set grid
set nokey

set multiplot

# 1ST
set size 0.5,0.5
set origin 0,0.5
set title "Btemp vs. Bperp" 
set xlabel "Btemp [days]"
set ylabel "Bperp [m]"
plot "INFILE" using 7:8

# 2ND
set size 0.5,0.5
set origin 0,0
set title "Bperp vs. Delta_fDC" 
set xlabel "Bperp [m]"
set ylabel "Delta_fDC [Hz]"
plot "INFILE" using 8:9

# 3RD
set size 0.5,0.5
set origin 0.5,0.5
set title "Btemp vs. Delta_fDC" 
set xlabel "Btemp [days]"
set ylabel "Delta_fDC [Hz]"
plot "INFILE" using 7:9

# 4TH -----------
set size 0.5,0.5
set origin 0.5,0
set title "Btemp vs. Bperp vs. Delta_fDC" 
set ticslevel 0.5                       # does this actually work?
set view 60,30
set xlabel "Btemp [days]"
set ylabel "Bperp [m]"
# set zlabel "Delta_fDC [Hz]"           # overlaps with plot 0.5 0.5
# zlabel manually set, since zlabel rotate doesnt work
set label 1 "Delta_fDC [Hz]" center rotate by 90 at graph 0, graph 0, graph 0.5 offset -7
splot "INFILE" using 7:8:9

unset multiplot
reset



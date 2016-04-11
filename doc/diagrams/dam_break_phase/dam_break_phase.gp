set terminal epslatex color
set output "dam_break_phase.tex"
set key top left

xmax = 3.6
set xrange [0:xmax]
set xzeroaxis
set xtics axis 0, 1.0, 3.0
set xlabel "$h$" offset 24,6,0

set yzeroaxis
set ytics axis 0, 1.0, 3.0
set ylabel "$hu$" offset 4,10,0

set arrow from 0,0 to 0,xmax size 0.1,15,60 \
  filled lc 0 lt 1
set arrow from 0,0 to xmax,0 size 0.1,15,60 \
  filled lc 0 lt 1

unset border

set style fill transparent solid 0.1 noborder

upper(x) = x*(2*sqrt(3) - 2*sqrt(x))
lower(x) = x*(-2*sqrt(3) + 2*sqrt(x))
plot '+' using 1:(upper($1)):(lower($1)) with filledcurves \
  above x2=3 notitle lc rgb "#00FF00", \
  '+' using 1:(upper($1)) title "$w_{1,\\max}$" lc 0 lt 1 with lines,\
  '+' using 1:(lower($1)) title "$w_{2,\\min}$" lc 0 lt 2 with lines,\
  "solution_exact.gpl" using 2:3 title "$\\mathbf{u}(x,t)$, Exact" lc 0 with lines,\
  "solution_initial.gpl" using 2:3 title "$\\mathbf{u}(x,0)$" lc 0 pt 4 with points,\
  "solution_Gal_FE.gpl" using 2:3 title "$\\mathbf{u}(x,t)$, Galerkin" lc 1 pt 1 with points,\
  "solution_DIV_FE.gpl" using 2:3 title "$\\mathbf{u}(x,t)$, DI" lc 2 pt 2 with points

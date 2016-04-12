set terminal epslatex color
set output "dam_break_phase.tex"

set key Left reverse top left

xmax = 3.6
set xrange [0:xmax]
set yrange [:3.3]

unset border

set xzeroaxis
set xtics axis 0, 1.0, 3.0
set xlabel "$h$" offset 23,6,0

set yzeroaxis
set ytics axis 0, 1.0, 3.0
set ylabel "$hu$" offset 3.5,11,0

set arrow from 0,0 to 0,xmax size 0.1,15,60 \
  filled lc 0 lt 1
set arrow from 0,0 to xmax,0 size 0.1,15,60 \
  filled lc 0 lt 1

set grid noxtics front

set style fill transparent solid 0.1 noborder

upper(x) = x*(2*sqrt(3) - 2*sqrt(x))
lower(x) = x*(-2*sqrt(3) + 2*sqrt(x))
plot '+' using 1:(upper($1)):(lower($1)) with filledcurves \
  above x2=3 notitle lc "blue", \
  '+' using 1:(upper($1)) title "$w_1(\\mathbf{u}) = \\max_i w_1(\\mathbf{u}^0_i)$" lc 0 lt 2 with lines,\
  '+' using 1:(lower($1)) title "$w_2(\\mathbf{u}) = \\min_i w_2(\\mathbf{u}^0_i)$" lc 0 lt 3 with lines,\
  "solution_exact.gpl" using 2:3 title "$\\mathbf{u}(x,t)$, Exact" \
    lc 0 lt 1 with lines,\
  "solution_initial.gpl" using 2:3 title "$\\mathbf{u}(x,0)$" \
    lc 4 lt 1 pt 6 ps 2 with linesp,\
  "solution_DIV_FE.gpl" using 2:3 title "$\\mathbf{u}(x,t)$, DI" \
    lc 2 lt 1 pt 2 with linesp,\
  "solution_Gal_FE.gpl" using 2:3 title "$\\mathbf{u}(x,t)$, Galerkin" \
    lc 1 lt 1 pt 1 with linesp

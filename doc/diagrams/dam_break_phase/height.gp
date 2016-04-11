set terminal epslatex color
set output "dam_break_height.tex"

set key Left reverse bottom left

set xrange [-5:5]

set ylabel "$h$" offset 2,0,0
set xlabel "$x$"
plot \
  "solution_exact.gpl" using 1:2 title "$\\mathbf{u}(x,t)$, Exact" \
    lc 0 lt 1 with lines,\
  "solution_Gal_FE.gpl" using 1:2 title "$\\mathbf{u}(x,t)$, Galerkin" \
    lc 1 lt 1 pt 1 with linesp,\
  "solution_DIV_FE.gpl" using 1:2 title "$\\mathbf{u}(x,t)$, DI" \
    lc 2 lt 1 pt 2 with linesp

set terminal epslatex color
set output "dam_break_phase.tex"
set ylabel "$hu$"
set xlabel "$h$"
set key top left

plot [0:4] x*(2*sqrt(3) - 2*sqrt(x)) title "$w_{1,\\max}$" lc 0 lt 1 with lines,\
  x*(-2*sqrt(3) + 2*sqrt(x)) title "$w_{2,\\min}$" lc 0 lt 2 with lines,\
  "solution_Gal_FE.gpl" using 2:3 title "$\\mathbf{u}(x)$, Galerkin" lc 1 pt 1 with points,\
  "solution_DIV_FE.gpl" using 2:3 title "$\\mathbf{u}(x)$, Domain-Invariant" lc 2 pt 2 with points

function [zq,wq] = get_GL_quadrature(nq)
% GLNodeWt  Nodes and weights for Gauss-Legendre quadrature of arbitrary order
%           obtained by solving an eigenvalue problem
%
% Synopsis:  [zq,wq] = GetGLQuadrature(nq)
%
% Input:     nq = number of quadrature points
%
% Output:    zq = vector of nodes
%            wq = vector of weights

%  Algorithm based on ideas from Golub and Welsch, and Gautschi.  For a
%  condensed presentation see H.R. Schwarz, "Numerical Analysis: A
%  Comprehensive Introduction," 1989, Wiley.  Original MATLAB
%  implementation by H.W. Wilson and L.H. Turcotte, "Advanced Mathematics
%  and Mechanics Applications Using MATLAB," 2nd ed., 1998, CRC Press

beta    = (1:nq-1)./sqrt(4*(1:nq-1).^2 - 1);
J       = diag(beta,-1) + diag(beta,1);    % eig(J) needs J in full storage
[V,D]   = eig(J);
[zq,iu] = sort(diag(D));  %  nodes are eigenvalues, which are on diagonal of D
wq      = 2*V(1,iu)'.^2;  %  V(1,ix)' is column vector of first row of sorted V

return
end
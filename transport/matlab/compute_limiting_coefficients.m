function L = compute_limiting_coefficients(F,uL,ML,W_max,W_min)
% computes the limiting coefficient matrix L
%
% L     = limiting coefficient matrix
% F     = flux correction matrix
% uL    = low-order solution
% W_max = upper bound for max principle
% W_min = lower bound for max principle

n=length(uL);

L      = zeros(n,n);
Pplus  = zeros(n,1);
Pminus = zeros(n,1);
Rminus = zeros(n,1);
Rplus  = zeros(n,1);
Qminus = zeros(n,1);
Qplus  = zeros(n,1);

% P vectors
for i=1:n   
    Pplus(i) =sum(max(0,F(i,:)));
    Pminus(i)=sum(min(0,F(i,:)));
end

% Q vectors
for i=1:n
    % changed to Zalesak definition (modified for steady-state)
    Qplus(i)  = ML(i,i)*(W_max(i) - uL(i));
    Qminus(i) = ML(i,i)*(W_min(i) - uL(i));
end

% R vectors
for i = 1:n
    if (Pplus(i) ~= 0)
        Rplus(i) = min(1, Qplus(i)/Pplus(i));
    else
        Rplus(i) = 1.0;
    end
    if (Pminus(i) ~= 0)
        Rminus(i) = min(1, Qminus(i)/Pminus(i));
    else
        Rminus(i) = 1.0;
    end
end
% set R+ and R- = 1 for Dirichlet node as Kuzmin recommended. This
% prevents L_ij from automatically being 0 for j in the support of i
Rplus(1)  = 1.0;
Rminus(1) = 1.0;

% limiting coefficients
for i=1:n
    for j=1:n
        if(F(i,j)>=0)
            L(i,j) = min(Rplus(i),Rminus(j));
        else
            L(i,j) = min(Rminus(i),Rplus(j));
        end
    end
end

return
end
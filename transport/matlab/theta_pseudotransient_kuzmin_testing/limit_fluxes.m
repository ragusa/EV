function Flim = limit_fluxes(F,u_aux,ML)
% computes the limiting coefficient matrix L
%
% L     = limiting coefficient matrix
%
% F     = flux correction matrix
% ML    = lumped mass matrix
% W_max = upper bound for max principle
% W_min = lower bound for max principle

n = length(u_aux);

Flim   = zeros(n,n);
Pplus  = zeros(n,1);
Pminus = zeros(n,1);
Rminus = zeros(n,1);
Rplus  = zeros(n,1);
Qminus = zeros(n,1);
Qplus  = zeros(n,1);

% P vectors
for i=1:n   
    Pplus(i) = sum(max(0,F(i,:)));
    Pminus(i)= sum(min(0,F(i,:)));
end

% Q vectors
for i = 1:n
    % compute index range of support of i
    i1=max(i-1,1);
    i2=min(i+1,n);
    % compute max and min old solution in the support of each dof
    u_max = max(u_aux(i1:i2));
    u_min = min(u_aux(i1:i2));
    % compute Q+/-
    Qplus(i)  = ML(i,i)*(u_max - u_aux(i));
    Qminus(i) = ML(i,i)*(u_min - u_aux(i));
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
Rplus(end) = 1.0;
Rminus(end)=1.0;

% limiting coefficients
for i=1:n
    for j=1:n
        % symmetric limitation
        if(F(i,j)>=0)
            Flim(i,j) = min(Rplus(i),Rminus(j))*F(i,j);
        else
            Flim(i,j) = min(Rminus(i),Rplus(j))*F(i,j);
        end
%         % upwind-biased limitation
%         k = min(i,j); % upwind node
%         if (F(i,j)>=0)
%             Flim(i,j) = Rplus(k)*F(i,j);
%         else
%             Flim(i,j) = Rminus(k)*F(i,j);
%         end
    end
end

return
end
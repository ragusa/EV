function LIM=fct(uH,D,uL,K)

% definition: f_ij = d_ij (u_i-u_j), j/=i

n=length(uH);
F=0*D;

for i=1:n
    for j=1:n
        F(i,j)=D(i,j)*(uH(i)-uH(j));
    end
end
check=abs(F+F');
if( max(max(check)) > 1e-10 )
    error('we should have F+F^t=0');
end

Pplus =zeros(n,1);
Pminus=zeros(n,1);
Rminus=zeros(n,1);
Rplus =zeros(n,1);
Qminus=zeros(n,1);
Qplus =zeros(n,1);
u_max =zeros(n,1);
u_min =zeros(n,1);

% compute lower/upper bounds
for i=1:n
    i1=max(i-1,1);
    i2=min(i+1,n);
    u_max(i) = max(uL(i1:i2));
    u_min(i) = min(uL(i1:i2));
end

% P vectors
for i=1:n
    for j = 1:n
        if (K(i,j) <= K(j,i))
            Pplus(i)  = Pplus(i)  + max(0,F(i,j));
            Pminus(i) = Pminus(i) + min(0,F(i,j));
        end
    end
    %Pplus(i) =sum(max(0,F(i,:)));
    %Pminus(i)=sum(min(0,F(i,:)));
end

% Q vectors
for i=1:n
    Qplus(i) = D(i,:)*max(0,uH-uH(i));
    Qminus(i)= D(i,:)*min(0,uH-uH(i));
end

% R vectors
% Rplus =min(1,Qplus ./Pplus );
% Rminus=min(1,Qminus./Pminus);
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

for i=1:n
    for j=1:n
        if(F(i,j)>=0)
            LIM(i,j)=Rplus(i);
        else
            LIM(i,j)=Rminus(i);
        end
    end
end

return
end
function LIM=fct(uH,MC,ML,D,uL)

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
    Pplus(i) =sum(max(0,F(i,:)));
    Pminus(i)=sum(min(0,F(i,:)));

end

% Q vectors
for i=1:n
    Qplus(i) = D(i,:)*max(0,uH-uH(i));
    Qminus(i)= D(i,:)*min(0,uH-uH(i));
end

% R vectors
Rplus =min(1,Qplus ./Pplus );
Rminus=min(1,Qminus./Pminus);

for i=1:n
    for j=1:n
        if(F(i,j)>0)
            LIM(i,j)=Rplus(i);
        else
            LIM(i,j)=Rplus(i);
        end
    end
end
% check if LIM is symmetric

return
end
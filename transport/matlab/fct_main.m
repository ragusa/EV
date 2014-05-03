close all; clc;

len=10;
nel=10;
omega=.1;
sigma=1;
src=0;
inc=1;
x=linspace(0,len,nel+1);

[MC,K,ML,D,b]=build_matrices(len,nel,omega,sigma,src,inc);

% low order solve
A=-(K+D);
A(1,:)=0;A(1,1)=1;
uL=A\b;
plot(x,uL,'s-'); hold all

% high order solve
A=-K;
A(1,:)=0;A(1,1)=1;
uH=A\b;
plot(x,uH,'+-')

% exact 
xx=linspace(0,len,1000);
if(sigma>eps)
    exact=src/sigma+(inc-src/sigma)*exp(-sigma*xx/omega);
else
    exact=inc+src/omega*xx;
end
plot(xx,exact,'r-')
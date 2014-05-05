clear all; close all; clc;

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
AL=A;
uL=AL\b;
plot(x,uL,'s-'); hold all
legend_entries=char('low-order');

% high order solve
A=-K;
A(1,:)=0;A(1,1)=1;
AH=A;
uH=AH\b;
plot(x,uH,'+-')
legend_entries=char(legend_entries,'high-order');

% compute limiter
LIM=fct(uH,MC,ML,D,uL);

% high order solve with limiter
A=-(K+D);
A(1,:)=0;A(1,1)=1;
LIM(:,:)=.25;
aux=-((LIM.*D)*uH);
aux(1)=0;
rhs=b+aux;
uLim=A\rhs;
plot(x,uLim,'v-')
legend_entries=char(legend_entries,'limited');

% exact 
xx=linspace(0,len,1000);
if(sigma>eps)
    exact=src/sigma+(inc-src/sigma)*exp(-sigma*xx/omega);
else
    exact=inc+src/omega*xx;
end
plot(xx,exact,'m-')
legend_entries=char(legend_entries,'exact');
legend(legend_entries)

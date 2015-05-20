function plot_viscosity(x,viscL,viscE)

x_cell = 0.5*(x(1:end-1)+x(2:end));
viscH = min(viscL,viscE);

semilogy(x_cell,viscL,'r-+');
hold on;
semilogy(x_cell,viscE,'b-x');
semilogy(x_cell,viscH,'g-o');
legend('Low','Entropy','High');
xlabel('x');
ylabel('Viscosity');
hold off;

end
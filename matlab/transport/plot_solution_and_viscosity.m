function plot_solution_and_viscosity(x,u,viscL,viscE)

figure

subplot(2,1,1);
plot(x,u,'k-x');
xlabel('x');
ylabel('Solution');
grid on;

subplot(2,1,2);
plot_viscosity(x,viscL,viscE);
grid on;

end
function plot_FCT(x,uH,uFCT,Wminus,Wplus)

figure

plot(x,Wminus,'k:o');
hold on;
plot(x,Wplus,'k:x');
plot(x,uH,'b-+');
plot(x,uFCT,'g-s');
hold off;
legend('W-','W+','High','FCT','Location','Best');

end
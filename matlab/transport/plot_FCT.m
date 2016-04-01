function plot_FCT(x,uH,uFCT,Wminus,Wplus,out_opts,new_figure,make_pause)

% unpack options
legend_location = out_opts.legend_location;
pause_type = out_opts.pause_type;
pausetime = out_opts.pausetime;

% if specified, create new figure
if (new_figure)
    figure(1);
    clf;
    hold on;
end

% plot sets
plot(x,Wminus,'b--');
plot(x,Wplus,'r--');
plot(x,uH,'b-+');
plot(x,uFCT,'g-s');

% legend
legend('W-','W+','High','FCT','Location',legend_location);

% pause if specified
if (make_pause)
    if (strcmp(pause_type,'wait'))
        k = waitforbuttonpress;
    else
        pause(pausetime);
    end
end

end
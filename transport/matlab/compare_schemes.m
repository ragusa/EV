function compare_schemes()
clear; clc; close all;

%% User options

plot_low_order  = false;
plot_high_order = false;
plot_FCT        = true;

include_Gal = true;
include_EV  = true;

include_FE      = true;
include_CN      = false;
include_BE      = false;
include_SSPRK33 = true;

legend_location = 'NorthEast';

save_plot = true;
output_file = 'output/stairstepping2.pdf';

%% Plot

% plot exact solution
alllegends = plot_set('output/uexact.txt','Exact','k-',[],[],[],true,false);

hold on;

% plot low-order solution
if (plot_low_order)
    if (include_FE)
        alllegends = plot_set('output/uL_FE.txt','Low-Order, FE','r-o','r',5,alllegends,false,true);
    end
    if (include_CN)
        alllegends = plot_set('output/uL_CN.txt','Low-Order, CN','r-s','r',5,alllegends,false,true);
    end
    if (include_BE)
        alllegends = plot_set('output/uL_BE.txt','Low-Order, BE','r-d','r',5,alllegends,false,true);
    end
    if (include_SSPRK33)
        alllegends = plot_set('output/uL_SSPRK33.txt','Low-Order, SSPRK33','r-^','r',5,alllegends,false,true);
    end
end

% plot high-order solution
if (plot_high_order)
    if (include_Gal)
        if (include_FE)
            alllegends = plot_set('output/uH_Gal_FE.txt','High-Order Galerkin, FE','b-o','k',7,alllegends,false,true);
        end
        if (include_CN)
            alllegends = plot_set('output/uH_Gal_CN.txt','High-Order Galerkin, CN','b-s','k',7,alllegends,false,true);
        end
        if (include_BE)
            alllegends = plot_set('output/uH_Gal_BE.txt','High-Order Galerkin, BE','b-d','k',7,alllegends,false,true);
        end
        if (include_SSPRK33)
            alllegends = plot_set('output/uH_Gal_SSPRK33.txt','High-Order Galerkin, SSPRK33','b-^','k',7,alllegends,false,true);
        end
    end
    if (include_EV)
        if (include_FE)
            alllegends = plot_set('output/uH_EV_FE.txt','High-Order EV, FE','b-o','b',5,alllegends,false,true);
        end
        if (include_CN)
            alllegends = plot_set('output/uH_EV_CN.txt','High-Order EV, CN','b-s','b',5,alllegends,false,true);
        end
        if (include_BE)
            alllegends = plot_set('output/uH_EV_BE.txt','High-Order EV, BE','b-d','b',5,alllegends,false,true);
        end
        if (include_SSPRK33)
            alllegends = plot_set('output/uH_EV_SSPRK33.txt','High-Order EV, SSPRK33','b-^','b',5,alllegends,false,true);
        end
    end
end

% plot FCT solution
if (plot_FCT)
    if (include_Gal)
        if (include_FE)
            alllegends = plot_set('output/uFCT_Gal_FE.txt','FCT Galerkin, FE','b-','k',7,alllegends,false,true);
        end
        if (include_CN)
            alllegends = plot_set('output/uFCT_Gal_CN.txt','FCT Galerkin, CN','g-s','k',7,alllegends,false,true);
        end
        if (include_BE)
            alllegends = plot_set('output/uFCT_Gal_BE.txt','FCT Galerkin, BE','g-d','k',7,alllegends,false,true);
        end
        if (include_SSPRK33)
            alllegends = plot_set('output/uFCT_Gal_SSPRK33.txt','FCT Galerkin, SSPRK33','c--x','k',7,alllegends,false,true);
        end
    end
    if (include_EV)
        if (include_FE)
            alllegends = plot_set('output/uFCT_EV_FE.txt','FCT EV, FE','r-','g',5,alllegends,false,true);
        end
        if (include_CN)
            alllegends = plot_set('output/uFCT_EV_CN.txt','FCT EV, CN','g-s','g',5,alllegends,false,true);
        end
        if (include_BE)
            alllegends = plot_set('output/uFCT_EV_BE.txt','FCT EV, BE','g-d','g',5,alllegends,false,true);
        end
        if (include_SSPRK33)
            alllegends = plot_set('output/uFCT_EV_SSPRK33.txt','FCT EV, SSPRK33','m--+','m',5,alllegends,false,true);
        end
    end
end

% add legend
legend(alllegends,'Location',legend_location);

% add annotations
xlabel('x');
ylabel('Solution');

% save plot
if (save_plot)
    export_fig(output_file,'-pdf','-transparent');
end

end
%% function to plot a single set
function alllegends = plot_set(filename,thislegend,format,marker_color,marker_size,alllegends,first,has_symbols)

% read data and plot
data = dlmread(filename);
x = data(:,1);
y = data(:,2);
if (has_symbols)
    plot(x,y,format,'MarkerFaceColor',marker_color,'MarkerEdgeColor',marker_color,'MarkerSize',marker_size);
else
    plot(x,y,format);
end

% add legend entry to list of legend entries
if (first)
    alllegends = char(thislegend);
else
    alllegends = char(alllegends,thislegend);
end

end
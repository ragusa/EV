function manual_limiting_ss_gui(x,uL,AL,inc,P,xexact,uexact,strong_dirichlet)

% number of antidiffusive fluxes
n = length(P);

% initialize limiting coefficients
L = zeros(n,1);
flim = zeros(length(uL),1);

% create figure
fig_x = 500;
fig_y = 300;
fig_height = 700;
fig_width = 900;
f = figure('Visible','off','Position',[fig_x,fig_y,fig_width,fig_height]);
movegui(f,'center');

% padding
padding = 10;
text_padding = 20;

% create table
col_names = {'Flux','Limiter'};
t = uitable(f,'Data',[P,L],'ColumnName',col_names,...
    'ColumnEditable',[false,true],'ColumnWidth',{75},...
    'CellEditCallback',@cell_edit,...
    'CellSelectionCallback',@cell_select);
table_width  = t.Extent(3);
table_height = min(t.Extent(4),fig_height-2*padding);
t.Position(3) = table_width;
t.Position(4) = table_height;
t.Position(1) = padding;
t.Position(2) = fig_height - table_height;

% create slider from 0 to 1
slider_width = 20;
slider_height = 100;
slider_x = t.Position(1) + table_width + padding;
slider_y = fig_height - slider_height;
myslider = uicontrol('Style','slider',...
    'Position',[slider_x,slider_y,slider_width,slider_height],...
    'String','Limiting Coefficient','Callback',@slider_change);

% create plot axes
axes_x = myslider.Position(1) + slider_width + padding + 2*text_padding;
axes_y1 = padding + 2*text_padding;
axes_width = fig_width - (table_width + slider_width + 3*padding + 4*text_padding);
axes_height = (fig_height - (3*padding + 4*text_padding))*0.5;
axes_y2 = axes_y1 + axes_height + padding + 2*text_padding;
axes1 = axes('Units','Pixels','Position',[axes_x,axes_y1,axes_width,axes_height]);
axes2 = axes('Units','Pixels','Position',[axes_x,axes_y2,axes_width,axes_height]);
hold(axes1,'on');
hold(axes2,'on');

% temporary handles
temp1 = [];
temp2 = [];
temp3 = [];
temp4 = [];

% create initial plot
uFCT = uL;
plot_data();

% make figure visible
set(f, 'Visible', 'on');

    % function to solve for FCT solution
    function compute_fct_solution()
        flim = zeros(length(uL),1);
        for i = 1:n
            flim(i) = flim(i) - P(i)*L(i);
            flim(i+1) = flim(i+1) + P(i)*L(i);
        end
        system_rhs = AL*uL + flim;
        system_matrix = AL; 
        if (strong_dirichlet)
            system_rhs(1) = inc;
            system_matrix(1,:) = 0; system_matrix(1,1) = 1;
        end
        uFCT = system_matrix \ system_rhs;
    end

    % function to plot data
    function plot_data
        cla(axes1);
        plot(axes1,xexact,uexact);
        plot(axes1,x,uL);
        plot(axes1,x,uFCT,'-x');
        xlabel(axes1,'x');
        ylabel(axes1,'Solution');
        legend(axes1,'Exact','Low','FCT');

        cla(axes2);
        plot(x,flim);
        xlabel('x');
        ylabel('Flux');
    end

    % callback function for selecting cell(s)
    function cell_select(src,evnt)
        % set user data
        set(src,'UserData',evnt.Indices);
        
        % get selection indices
        ind = get(t,'UserData');
        
        % put arrow symbols on selected nodes
        n_selections = length(ind(:,1));
        if (n_selections == 1)
            row = ind(1,1);
            delete(temp1);
            delete(temp2);
            delete(temp3);
            delete(temp4);
            if (P(row) > 0)
                temp1 = plot(axes2,x(row),flim(row),'bv');
                temp2 = plot(axes2,x(row+1),flim(row+1),'r^');
                temp3 = plot(axes1,x(row),uFCT(row),'bv');
                temp4 = plot(axes1,x(row+1),uFCT(row+1),'r^');
            else
                temp1 = plot(axes2,x(row),flim(row),'r^');
                temp2 = plot(axes2,x(row+1),flim(row+1),'bv');
                temp3 = plot(axes1,x(row),uFCT(row),'r^');
                temp4 = plot(axes1,x(row+1),uFCT(row+1),'bv');
            end
        end
    end

    % callback function for editing a table cell
    function cell_edit(~,~)
        % update L
        data = get(t,'Data');
        L = data(:,2);
        
        % recompute FCT solution
        compute_fct_solution();
        
        % replot
        plot_data();
    end

    % callback function for changing slider
    function slider_change(~,~)
        % get selection indices
        ind = get(t,'UserData');
        
        % if at least one row has been selected
        if (~isempty(ind))
            % get slider value
            value = get(myslider, 'Value');
            
            % loop over rows of selection
            n_selections = length(ind(:,1));
            for i = 1:n_selections
                % get row index
                row = ind(i,1);
                
                % change value in data array
                L(row) = value;
            end
                                    
            % update table
            set(t,'Data',[P,L]);
            
            % recompute FCT solution
            compute_fct_solution();
            
            % replot
            plot_data();
        end
    end

end

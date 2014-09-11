close all; clc; clear;

% CFLcases = char('0.95','1.5','3.0','5.0','10.0');
CFLcases = char('0.1','0.25','0.5','0.75','0.95');
imex = 'ex';

ncases = size(CFLcases,1);

colorlist = ['k','r','g','b','m'];
figure('units','normalized','position',[.1 .1 .4 .4]);
xlabel('h');
ylabel('Discrete L_2 Error');
set(gca,'xscale','log');
set(gca,'yscale','log');
hold on;
for i = 1:ncases
    CFLstring = deblank(CFLcases(i,:));
    filename = ['convergence_',imex,CFLstring,'.csv'];
    convdata = csvread(filename);
    h      = convdata(:,1);
    errL   = convdata(:,2);
    errH   = convdata(:,3);
    errFCT = convdata(:,4);
    formatL   = [colorlist(i),'--+'];
    formatH   = [colorlist(i),'--o'];
    formatFCT = [colorlist(i),'--x'];
    loglog(h,errL,formatL);
    loglog(h,errH,formatH);
    loglog(h,errFCT,formatFCT);
    legendL   = ['CFL = ',CFLstring,', Low'];
    legendH   = ['CFL = ',CFLstring,', High'];
    legendFCT = ['CFL = ',CFLstring,', FCT'];
    if (i==1)
        legendstrings = char(legendL,legendH,legendFCT);
    else
        legendstrings = char(legendstrings,legendL,legendH,legendFCT);
    end
end

href = [1,10/(10*2^4)];
e1 = errL(1);
errref = [e1,e1*(href(2)/href(1))];
loglog(href,errref,'c-');
legendstrings = char(legendstrings,'m = 1');
e1 = errH(1);
errref = [e1,e1*(href(2)/href(1))^4];
loglog(href,errref,'c-');
legendstrings = char(legendstrings,'m = 4');

legend(legendstrings,'Location','EastOutside');

export_fig('explicit_lowCFL.pdf','-transparent','-pdf');
% export_fig('implicit_lowCFL.pdf','-transparent','-pdf');
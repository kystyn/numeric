function main()
    format long;
    rng('shuffle');
    [conds, ev] = prepareSLAEPlotCond();
    system('lab2.exe');
    disp('C++ done');
    drawPlot(conds, 'Dependence of x error on matrix cond, b error = const', 'cond_{M}');
    
    %[conds ev] = prepareSLAEPlotXBerr();
    %system('lab2.exe');
    %disp('C++ done');
    %drawPlot(conds * ev, 'Dependence of x error on b error, cond M = const');
end
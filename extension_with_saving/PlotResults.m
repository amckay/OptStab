function PlotResults()


FM = load('/Users/iragm01/projects/OptStab/WithSaving/Code/FullModel/results/WorkTable.mat');
NS = load('/Users/iragm01/projects/OptStab/WithSaving/Code/NoSavings/results/WorkTable.mat');


FM.WorkTable(:,11) = FM.WorkTable(:,11) * 100;
NS.WorkTable(:,11) = NS.WorkTable(:,11) * 100;


f2 = InnerPlot(11)
ylabel('St. dev. of unemployment (\times 10^{-2})')



function f = InnerPlot(itoplot)


[~,I] = sort(FM.WorkTable(:,7));
FM.WorkTable = FM.WorkTable(I,:);

[~,I] = sort(NS.WorkTable(:,7));
NS.WorkTable = NS.WorkTable(I,:);


f = figure();
p = plot(NS.WorkTable(:,7),NS.WorkTable(:,itoplot),FM.WorkTable(:,7),FM.WorkTable(:,itoplot));
lg = legend('Baseline', 'Extended Model With Savings','Location','NorthEast');

sizes = {2,2,1};
styles = {'--','-','-'};
markers = {'none','none','+'};
for i = 1:2
    p(i).LineWidth = sizes{i};
    p(i).LineStyle = styles{i};
    p(i).Color = 'k';
    p(i).Marker = markers{i};
end

lg.FontSize = 14;
lg.Box = false;
f.CurrentAxes.FontSize = 14;
xlim([0.063, 0.068]);

xlabel('Steady state unemployment rate')
% yticks([0.014:0.001:0.017])
end

end
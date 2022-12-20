% Figure S13

parcs = [1 5 3 6]; 

figure('DefaultAxesFontSize', 14); 
t = tiledlayout(2, 2);
for ii = 1:4
    nexttile();
    histogram(mean(PARC_SA{parcs(ii)}), 20); 
    title(parc_name2(parcs(ii)));
end
xlabel(t, 'Node Surface Area (mm^2)', 'FontSize', 18);
ylabel(t, 'Count', 'FontSize', 18);


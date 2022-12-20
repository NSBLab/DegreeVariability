function sa = viewParcellation(verts, faces, parc)

figure;%("Position", [2000 600 1200 600]); 
tiledlayout(2, 3);

nexttile(1);
bar(ROI_Nverts(parc));
ylabel("N verts");

nexttile(4);
sa = calcRoiArea(verts, faces, parc);
bar(sa);
ylabel("SA");
xlabel("ROI ID");

nexttile(2, [2 2]);
brainplot(verts, faces, parc, sa, 0, parula(100), 'SurfaceArea');

end
PARCNAMES = {'aparc','random200','HCPMMP1','random500','Schaefer200_17net','Schaefer500_17net'};
addpath /usr/local/freesurfer/5.3/matlab/
for p = 1:6
PARCNAME = PARCNAMES{p};
clear ROISA
for i = 1:length(subjects)

    SUBID = subjects(i);
    
    for j = 1:2
        if j == 1
            hemi = 'l';
        else
           hemi = 'r'; 
        end
    
parc_path = ['/fs03/kg98/stuarto/GenCog_xnat/MRH035_',num2str(SUBID),'_MR01/label/',hemi,'h.',PARCNAME,'.annot'];

surf_path = ['/fs03/kg98/stuarto/GenCog_xnat/MRH035_',num2str(SUBID),'_MR01/surf/',hemi,'h.white'];

[vertices, label, colortable] = read_annotation(parc_path);

VERT_FSids = colortable.table(:,5);

Nrois = length(VERT_FSids);

P = changem(label,0:Nrois-1,VERT_FSids);
if p == 3
    P(P==16777215) = 0;
end
[V, F] = read_surf(surf_path);

ROISA_temp{j} = GetParcSurfArea(F+1,V,P);
    end
    
ROISA(i,:) = [ROISA_temp{1}(2:end); ROISA_temp{2}(2:end)];

end

PARC_SA{p} = ROISA;

end
function [ROISA,ROI] = GetParcSurfArea(F,V,P)

% Calculate surface area
a = V(F(:, 2), :) - V(F(:, 1), :);
b = V(F(:, 3), :) - V(F(:, 1), :);
c = cross(a, b, 2);
areaTri = 1/2 *(sqrt(sum(c.^2, 2)));

% Get ROI id for each vertex in each face
Fp = P(F);

% Find faces that exist within a ROI
WithinROI = max(Fp,[],2)-min(Fp,[],2)==0;

% Assign each face to a ROI, ignoring those which have split membership
F_ROI_ID = Fp(:,1);
F_ROI_ID(~WithinROI)= NaN;

% Get the unique list of ROIs
ROI = unique(F_ROI_ID(WithinROI))';

ROISA = zeros(length(ROI),1);

for i = 1:length(ROI)
   
    ROI_F_id = F_ROI_ID == ROI(i);
    ROISA(i) = sum(areaTri(ROI_F_id));
    
end
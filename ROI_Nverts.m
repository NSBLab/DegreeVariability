function out = ROI_Nverts(parc, includeZero)
if nargin < 2
    includeZero = false;
end

%lowerlimit  = 0.5 - includeZero;
%lowerlimit = tern(includeZero, -0.5, 0.5);
if includeZero
    lowerlimit = -0.5;
else
    lowerlimit = 0.5;
end

out = histcounts(parc, lowerlimit:(max(parc)+0.5));
% out = histcounts(parc, tern(includeZero, -0.5, 0.5):(max(parc)+0.5));

end

% DESCRIPTION: Outputs the global energy minima

function outData = ComputeMins(DVectors,IntCurves)

%Extract the interaction energy curves and its stagger vector: 
outData = zeros(size(DVectors,1),2);
OrderK = 1; %order of the minimum (i.e. 1 corresponds to the global minimum)

for idx = 1 : size(DVectors,1)
    %Extract the interaction curve data:
    ECurve = IntCurves{idx};
    DVec = DVectors{idx};
    [~,locs] = findpeaks(-ECurve,'SortStr','descend','NPeaks',OrderK);
    
    %Save the mean energy and associated std value: 
    outData(idx,:) = [ECurve(locs(OrderK)),DVec(locs(OrderK))];
end

end
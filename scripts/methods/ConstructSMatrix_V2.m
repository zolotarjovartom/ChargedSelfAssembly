% DESCRIPTION: Performs the same functionality as GetRobustMinima.m code
% without incorporation of robustness.

% NOTES: 

% SmoothWindow is a 2 element vector which corresponds to the smoothing
% interval (in Angstroms) around each discrete stagger value. 

function SMats = ConstructSMatrix_V2(DVectors_S,IntCurves_S,DVectors,IntCurves,...
    D,SpiralPairs,n,ns,SmoothWindow)
%Start by computing the minima of the smoothed curves: 

%Compute the first 3 minima for each interaction energy curve: 
NumSpiralPairs = numel(DVectors_S);
NumK = 3; %Number of minima

%Define the container for the minima: 
AllMins = zeros(NumSpiralPairs*NumK,1);

for idx = 1 : NumSpiralPairs
    if numel(IntCurves_S{idx}) > 3
        %Compute the minima: 
        [~,Locs] = findpeaks(-IntCurves_S{idx},'SortStr','Descend','NPeaks',NumK);
        
        %Identify the corresponding minimal energy in the unsmoothed potential:
      
        %Construct the search intervals: 
        SearchInts = [DVectors_S{idx}(Locs)'-SmoothWindow(1),...
            DVectors_S{idx}(Locs)'+SmoothWindow(2)];
        minCont = zeros(NumK,1);
        
        for minIdx = 1 : NumK
            %Identify the smoothing region: 
            SmoothReg = find(and(DVectors{idx} > SearchInts(minIdx,1),...
                DVectors{idx} < SearchInts(minIdx,2)));
    
            %Compute the minimal value of unsmoothed energy over SmoothReg:
            [~,minLoc] = min(IntCurves{idx}(SmoothReg));
    
            %Extract the corresponding stagger: 
            minCont(minIdx) = SmoothReg(minLoc);      
        end
    
        %Store the staggers:
        AllMins((3*idx - 2): 3*idx)  = DVectors{idx}(minCont)';

    else
        AllMins((3*idx - 2): 3*idx) = zeros(numel((3*idx - 2): 3*idx),1);
    end
end

%Populate the stagger matrix: 
SMats = GetAmpMatrix(AllMins,D,SpiralPairs,n,ns);
end

%SUBROUTINES: 
%Construct and save amplitude matrices based on the 
function [SMats] = GetAmpMatrix(Staggers,D,SpiralPairs,n,ns)
%Compute the stagger matrices:
SMats = PopulateMatrix(Staggers,SpiralPairs,1,ns); %Antisymmetric matrix

%For each matrix, calculate it modulo n*D: 
SMats = cellfun(@(X) mod(X,n*D), SMats,'UniformOutput',0);

%Define a subroutine to populate a matrix S_ij, where ij refer to the
%interacting pairs of strips:
function [PopMats] = PopulateMatrix(InArr,SpiralPairs,IsAntisym,ns)
    %Generate stagger matrices:
    SG = zeros(ns,ns); %Global minima
    SS = zeros(ns,ns); %Secondary minima
    ST = zeros(ns,ns); %Tertiary minima

    %Identify the global, secondary and tertiary minima: 
    Idx = sub2ind(size(SG),SpiralPairs(:,1),SpiralPairs(:,2)); %Convert into a linear index
    IdxRev = sub2ind(size(SG),SpiralPairs(:,2),SpiralPairs(:,1));
  
    %Populate the stagger matrices:
    SG(Idx) = InArr(1:3:end); SS(Idx) = InArr(2:3:end); ST(Idx) = InArr(3:3:end);
    
    if IsAntisym == 1 
        %Ensure that the matrix is antisymmetric:
        SG(IdxRev) = -1.*InArr(1:3:end); SS(IdxRev) = -1.*InArr(2:3:end); ST(IdxRev) = -1.*InArr(3:3:end);
        
        %Correct the diagonal elements (they do not need a minus sign):
        s = size(SG); index = 1:s(1)+1:s(1)*s(2);
        SG(index) = SG(index).*-1; SS(index) = SS(index).*-1; ST(index) = ST(index).*-1;
         
    elseif IsAntisym == 0 
        %Populate below the diagonal:
        SG(IdxRev) = InArr(1:3:end); SS(IdxRev) = InArr(2:3:end); ST(IdxRev) = InArr(3:3:end);
    
    end
    %Output the stagger matrices: 
    PopMats = {SG,SS,ST};
end
end
% DESCRIPTION: Outputs a matrix cell, wherein I(i,j) = E_{i-j}

function IMat = GetIMat(IntCurves,SVec,n,D,NumSpirals,SpiralPairs)
%Define containers:
ModIntData = cell(numel(IntCurves),1);

%Modify the interaction curve data: 
for cIdx = 1 : numel(IntCurves)    
   %Convert the energy to be on the interval [0,n*D):
   [SVec_C,uIdx] =  unique(mod(SVec{cIdx},n*D));
   
   %Extract the energy corresponding to the interaction with no noise:
   ModIntData{cIdx} = [IntCurves{cIdx}(1,uIdx);SVec_C; zeros(size(SVec_C))];

   %Append an entry at n*D: 
   ModIntData{cIdx} = horzcat(ModIntData{cIdx},...
              [ModIntData{cIdx}(1,1); n*D; 0]);
end

%Allocate the interaction curve data into an interaction matrix: 
IMat = cell(NumSpirals,NumSpirals);

for mIdx = 1 : numel(IntCurves)
    %Populate the diagonal and above diagonal entries:
    IMat{SpiralPairs(mIdx,1),SpiralPairs(mIdx,2)} = ModIntData{mIdx};
    
    %Populate the lower diagonal:
    if SpiralPairs(mIdx,1) ~= SpiralPairs(mIdx,2)
        %Rescale the stagger vector by -1: 
        [ModIntData_Sort,SortIdx] = unique(mod(-1.*ModIntData{mIdx}(2,:),n*D));
        ModIntData_Sort = vertcat(ModIntData{mIdx}(1,SortIdx),ModIntData_Sort,ModIntData{mIdx}(3,SortIdx));
       
        %Append an entry at n*D:
        ModIntData_Sort = horzcat(ModIntData_Sort,...
              [ModIntData_Sort(1,1); n*D; 0]);

        %Populate the lower diagonal entries:
        IMat{SpiralPairs(mIdx,2),SpiralPairs(mIdx,1)} =  ModIntData_Sort; 
    end
end

end

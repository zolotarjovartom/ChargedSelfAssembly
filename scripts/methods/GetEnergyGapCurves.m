% DESCRIPTION: Computes the microfibrillar energy and eqbm probabilities 
% as a function of the gap lengthscale 
function outStruc = GetEnergyGapCurves(pKaData,pH)
%Define the range of gap values:
MidPt = 450; %Rougly the size of gap region
gapRange = MidPt.*[0.8,1.2]; %In Angstroms

%Convert into the range of D values: 
L = 2883.8; %In Angstroms
n = 5; v = 5;
DRange = (gapRange+L)./n; step = 0.2;
DRange = min(DRange):step:max(DRange);

%Slice the energy as a function of pH data: 
E = cell(numel(pH),1);
for idx = 1 : numel(pH)
    E{idx} = cellfun(@(X) X(idx,:), pKaData.Energy,'UniformOutput',0);
end

%Obtain the staggers:
Stagger = pKaData.Staggers;

%Construct the output containers: 
pOut = cell(numel(pH),1); EqStates = cell(numel(pH),1);

%Loop over the different pH values: one possible place to parallelise over
for pidx = 1 : numel(pH)
   %Compute the energy/probability curves as a function of gap size: 
   [pOut{pidx},EqStates{pidx}] = ComputeGlobalMinimaV2(E{pidx},...
       Stagger,DRange,n,v);
end

%Define an output structure:
outStruc = struct('Probs',{pOut},'EqStates',{EqStates},'DRange',{DRange});
end
%DESCRIPTION: Computes interaction energy curves as a function of pH and
%varying pKa values. Conventionally, this script varies the values of the
%pKa's of Asp and Lys

function pKaGrid(SpeciesNames,Species,SVecs,StripAACoors,StripAAInts,ContactMat,ColName,FPathpKa)
%Fix the ionic parameters:
I = 150; 
a0 = 4e-10; %Distance of closest approach 
eps_r = 76.5; %relative permittivity, take value for water at 303.15K
T = 303.15; %Temperature in K 

%Fix the parameters for the construction of ionisation windows: 
sf = 1; %Intervals are disjoint
alphaBounds = [0.05,0.95];
delta_alpha = 0.05; %Maximum ionisation increment for a pH increment
stepSize = 0.05; %This is the step size for the non-overlap measure. The step size is the
%fraction of the ionisation window that is traversed per incement in the
%pKa value.

%Construct the modified pKa values: 
pKaMat = GetModified_pKa(sf,alphaBounds,I,stepSize,eps_r,T);

%Compute the overlap measure: (sufficient here to do for just one of the
%residue pairs, here Arg and Lys) and ionisation intervals 
[mu,WindowsArg,WindowsLys] = ComputeOverlap(1,4,1,1,"B","B",pKaMat,alphaBounds,eps_r,T,I);
[~,WindowsGlu,WindowsAsp] = ComputeOverlap(3,2,-1,-1,"A","A",pKaMat,alphaBounds,eps_r,T,I);

%Convert the ionisation intervals into pH ranges to discretise: 
[basicMin_pH,basicMax_pH] = bounds([WindowsArg,WindowsLys],2);
[acidicMin_pH,acidicMax_pH] = bounds([WindowsGlu,WindowsAsp],2);

%Compute the discretised pH intervals: 
ResType = {'B','A','A','B','A','A','B'}';
qvec = [1, -1, -1, 1, -1, -1, 1]';
basic_pHVec = arrayfun(@(X,Y,Z) discretise_pH([X,Y],qvec,I,eps_r,T,pKaMat(Z,:),ResType,...
    delta_alpha,alphaBounds),basicMin_pH,basicMax_pH,transpose(1:size(pKaMat,1)),...
    'UniformOutput',0);
acidic_pHVec = arrayfun(@(X,Y,Z) discretise_pH([X,Y],qvec,I,eps_r,T,pKaMat(Z,:),ResType,...
    delta_alpha,alphaBounds),acidicMin_pH,acidicMax_pH,transpose(1:size(pKaMat,1)),...
    'UniformOutput',0);

%Merge the pH intervals: 
pH_vec = cellfun(@(X,Y) [X',Y'],acidic_pHVec,basic_pHVec,'UniformOutput',0);

%Set other parameters for the calculation: 
NumSpirals = 7; SpiralSet = 1; 

%Extract the relevant species data: 
SpecLoc = strcmp(SpeciesNames,Species);
StripAACoors = StripAACoors(SpecLoc); StripAAInts = StripAAInts(SpecLoc);
SVecs = SVecs.(Species); %Staggers 

%Construct containers for the outputs: 
tempStagger = cell(size(pKaMat,1),1);
tempE = cell(size(pKaMat,1),1);

%Run a parallel loop across the different values of pKa:
%NumWorkers = 15; 
%parpool('Processes',NumWorkers)
%Change the loop below to be parallel.
for paridx = 1 : size(pKaMat,1)
    disp('Starting Calculation number:')
    disp(paridx)
    [tempStagger{paridx},tempE{paridx}] = GetInterCurves_Bulk_DH_NonLocal(...
        NumSpirals,SpiralSet,ContactMat,StripAACoors{1},StripAAInts{1},...
    SVecs,I,a0,pH_vec{paridx},pKaMat(paridx,:),eps_r,T);
    disp('Finished Calculation number:')
    disp(paridx)
end

%Construct output structures: 
RawNames = strcat('muMeasure_',strtrim(replace(cellstr(num2str(round(mu,3))),{'-','.'},...
    {'minus','_'})));
data = zeros(1,size(RawNames,1)); 
S = cell2struct(num2cell(data,1), RawNames, 2);
S.('pHVec') = pH_vec; S.('I') = I; S.('a0') = a0; S.('Species') = Species;
S.('muVec') = mu; S.('pKaMat') = pKaMat;
subStruc = struct('Energy',[],'Staggers','');

%Populate the output structure: 
for idx = 1 : size(pKaMat,1)
    %Populate the substructure:
    subStruc.Energy = tempE{idx}; subStruc.Staggers = tempStagger{idx};
    
    %Populate the main structure: 
    S.(RawNames{idx}) = subStruc;
end

%Save the output: 
FName = strcat(FPathpKa,Species,'|Type',ColName{1},'|',ColName{2},'|','pKaGrid.mat');
save(FName,'S','-v7.3');
end

%SUBROUTINES: 

%Compute the overlap measure for the given set of pKa's 

%The convention is that residue B is always the one whose pKa is being
%modified

function [mu,WindowsA,WindowsB] = ComputeOverlap(resA,resB,qA,qB,resTypeA,resTypeB,pKaMat,alphaBounds,...
    eps_r,T,I)
    %Start by computing the ionisation intervals:
    numRows = size(pKaMat,1);
    [WindowsA,~] = get_pHInterval(repmat(alphaBounds,numRows,1),...
        repmat(resTypeA,numRows,1),repmat(qA,numRows,1),eps_r,T,I,pKaMat(:,resA));
    [WindowsB,~] = get_pHInterval(repmat(alphaBounds,numRows,1),repmat(resTypeB,...
        numRows,1),repmat(qB,numRows,1),eps_r,T,I,pKaMat(:,resB));
    midPt = @(X) X(:,1) + 0.5.*diff(X,1,2); %Identify midpoint of the ionisation window
    
    %Compute the overlap:
    mu = OverlapMeasure(midPt(WindowsA),midPt(WindowsB),WindowsA,WindowsB);
end

%Outputs a matrix containing modified pKa values. Each row corresponds to a
%different degree of interval overlap and each column to one of the 7
%ionisable residues 

%CONVENTIONALLY, we modify the pKa of columns 2 and 4
function pKaMat = GetModified_pKa(sf,alphaBounds,I,stepSize,eps_r,T)
    %Construct non-overlapping ionisation windows:
    [pKa_ref,pKa_shift] = Modify_pKa(1,sf,alphaBounds,I,stepSize,eps_r,T);
    
    %Constuct a range of pKa values to test our model on: 
    
    %Convention: modify the 2nd and 4th pKa values across the calculations
    pKaMat1 = repmat(pKa_ref,numel(pKa_shift),1) - ...
        [zeros(numel(pKa_shift),1),sort(pKa_shift,'ascend')',zeros(numel(pKa_shift),1),...
        sort(pKa_shift,'ascend')',zeros(numel(pKa_shift),3)];
    
    %Now for the other ordering: 
    [pKa_ref,~] = Modify_pKa(2,sf,alphaBounds,I,stepSize,eps_r,T);
    pKaMat2 = repmat(pKa_ref,numel(pKa_shift),1) + ...
        [zeros(numel(pKa_shift),1),sort(pKa_shift,'descend')',zeros(numel(pKa_shift),1),...
        sort(pKa_shift,'descend')',zeros(numel(pKa_shift),3)];
    
    %Combine together the matrices: 
    pKaMat1(end,:) = []; %Remove full overlap as this will be duplicated otherwise
    pKaMat = vertcat(pKaMat1,pKaMat2);
end
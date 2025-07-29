% DESCRIPTION: Computes the equilibrium statistics for 5-membered collagen
% microfibrils at different values of pH and pKa. 

function runProbsOnpKaGrid
%Specify the path where you will save the output file:
FPathEqbm = "/Users/artomzolotarjov/Documents/ChargedGitHub/scripts/eqbmStatistics/" ;

%Add the "methods" scripts: 
addpath("/Users/artomzolotarjov/Documents/ChargedGitHub/scripts/methods/");

%Load the pKa data:
ColType = "II"; %Specify the collagen type: 
TrimerType = "A1|A1|A1"; %Specify the trimer type:
Species = 'BosTaurus';
FPathpKa = '/Volumes/CollagenDrive/Debye-Huckel Interaction Curves/Updated Data 09.23.2024/';
FName = strcat(FPathpKa,Species,'|Type',ColType,'|',TrimerType,'|','pKaGrid.mat');
pKaData = load(FName);

%Identify the fields with interaction data:
fields = fieldnames(pKaData.S);
TFData = contains(fields,'muMeasure'); 

%Slice the interaction data: 
intDataSliced = cell(sum(TFData),1); slicedFields = fields(TFData);
for idx = 1 : sum(TFData)
    intDataSliced{idx} = pKaData.S.(slicedFields{idx});
end

%Extract the pH data:
pHData = pKaData.S.pHVec;

%Define an output container: 
outCell = cell(sum(TFData),1);

%PARALLELISE HERE:
%Run a parallel loop on the different values of mu measure: 
%NumWorkers = 15;
%parpool('Processes',NumWorkers)
for paridx = 1 : sum(TFData) 
    disp('Starting Calculation number:')
    disp(paridx)
    outCell{paridx} = GetEnergyGapCurves(intDataSliced{paridx},pHData{paridx});
    disp('Finished Calculation number:')
    disp(paridx)
end

%Define an output structure: 
data = zeros(1,size(slicedFields,1)); 
S = cell2struct(num2cell(data,1), slicedFields, 2);

%Populate the output structure: 
for idx = 1 : sum(TFData)
    S.(slicedFields{idx}) = outCell{idx};
end

%Append the relevant parameters to the output structure: 
for idx = sum(TFData) + 1 : numel(fields)
    S.(fields{idx}) = pKaData.S.(fields{idx});
end

%Save the output: 
FName = strcat(FPathEqbm,'EqbmProbspKaGrid_','_',Species,'|Type',ColType,'|',TrimerType,'.mat');
save(FName,'S','-v7.3');
end
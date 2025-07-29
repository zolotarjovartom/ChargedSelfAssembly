% DESCRIPTION: Produces Fig.5 of "The Role of Ionic Conditions and Charged Residues in
% the Emergence of D-Banded Collagen Microfibrils"

function genMainFig5 
%Specify the path to the eqbm statistics data: 
FPathEqbm = '/Volumes/CollagenDrive/Debye-Huckel Interaction Curves/Updated Data 09.23.2024/EqbmProbspKaGrid_BosTaurus|TypeII|A1|A1|A1.mat';

%Add the "methods" directory to the path: 
addpath("/Users/artomzolotarjov/Documents/ChargedGitHub/scripts/methods/");

%Specify the path to the residue data: 
FPathRes = '/Volumes/CollagenDrive/Bulk Residue Data/BulkResData_TypeII|A1|A1|A1_Mammals.mat';

%Obtain the residue counts: 
speciesSplit = 'Bos taurus';
NumStrips = 7; %number of helical strips 
getResCounts(FPathRes,speciesSplit,NumStrips);

%Load the eqbm data: 
eqbmData = load(FPathEqbm);

%Run the calculation: 
CollagenPhaseDiagram(eqbmData);
end
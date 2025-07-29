% DESCRIPTION: Produces Fig.8 of "The Role of Ionic Conditions and Charged Residues in
% the Emergence of D-Banded Collagen Microfibrils"

function genMainFig8
%Specify the path to the eqbm statistics data: 
FPathEqbm = '/Volumes/CollagenDrive/Debye-Huckel Interaction Curves/Updated Data 09.23.2024/EqbmProbspKaGrid_BosTaurus|TypeII|A1|A1|A1.mat';

%Add the "methods" directory to the path: 
addpath("/Users/artomzolotarjov/Documents/ChargedGitHub/scripts/methods/");

%Load the eqbm data: 
eqbmData = load(FPathEqbm);

%Run the calculation:
PaxConservedFeaturesV3(eqbmData);
end
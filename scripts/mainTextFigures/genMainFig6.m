% DESCRIPTION: Produces Fig.6 of The "Role of Ionic Conditions and Charged Residues in
% the Emergence of D-Banded Collagen Microfibrils"

function genMainFig6
%Specify the path to the eqbm statistics data: 
FPathEqbm = '/Volumes/CollagenDrive/Debye-Huckel Interaction Curves/Updated Data 09.23.2024/EqbmProbspKaGrid_BosTaurus|TypeII|A1|A1|A1.mat';

%Add the "methods" directory to the path: 
addpath("/Users/artomzolotarjov/Documents/ChargedGitHub/scripts/methods/");

%Load the eqbm data: 
eqbmData = load(FPathEqbm);

%Run the script for acidic phase transitions at mu_ref=0.15:
PlotEqProb_EqEnergy_V2(eqbmData,24)

%Run the script at mu_ref=0.75: 
PlotEqProb_EqEnergy_V2(eqbmData,36)
end
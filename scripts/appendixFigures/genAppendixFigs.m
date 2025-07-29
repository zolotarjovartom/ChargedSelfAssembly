% DESCRIPTION: Produces appendix figures of "The Role of Ionic Conditions and Charged Residues in
% the Emergence of D-Banded Collagen Microfibrils".

function genAppendixFigs
%Add the "methods" directory to the path: 
addpath("/Users/artomzolotarjov/Documents/ChargedGitHub/scripts/methods/");

%Figure 10:
%{
%Specify the relevant paths: 
FPathRes = '/Volumes/CollagenDrive/Bulk Residue Data/';

%Specify the species for which to plot residue frequencies: 
speciesSplit = 'Bos taurus';

%Run the script:
chargedStats(FPathRes,speciesSplit)
%}

%Figure 11: 
%{
%Specify the path to the eqbm statistics data: 
FPathEqbm = '/Volumes/CollagenDrive/Debye-Huckel Interaction Curves/Updated Data 09.23.2024/EqbmProbspKaGrid_BosTaurus|TypeII|A1|A1|A1.mat';

%Load the data: 
eqbmData = load(FPathEqbm);

%Run the calculation: 
PlotEqStateDBand(eqbmData)
%}

end

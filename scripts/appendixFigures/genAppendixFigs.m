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

%Figure 12: 
%{
%Specify the path to the interaction energy curves: 
FPathE = '/Volumes/CollagenDrive/Debye-Huckel Interaction Curves/Updated Data 09.23.2024/BosTaurus|TypeII|A1|A1|A1|pKaGrid.mat';
Data = load(FPathE);

%Extract the relevant data: 
tValmu = 0.15; %Make sure this matches with the field when extracting E and S!
muVec = Data.S.muVec;
[~,minLoc] = min(abs(muVec-tValmu));
EMat = Data.S.muMeasure_0_15.Energy; SVecs = Data.S.muMeasure_0_15.Staggers;
pHVals = Data.S.pHVec; pHVec = pHVals{minLoc};

%Choose pH for which the acidic and basic residues are fully-ionised: 
pHLoc = 24;
disp(strcat('Selected pH value is: ',num2str(pHVec(pHLoc))));

%Extract the relevant interaction energy curves: 
EVecs = cellfun(@(X) X(pHLoc,:), EMat,'UniformOutput',0);

%Compute the periodic interaction potential minima:
n = 5; 
D = (335.8 + 3000)/n; 

%Define parameters for the calculation of periodic interaction energy curves:
Window = [14,14]; %Smoothing is only used for finding minima  
SSpacing = 0.5;

%Compute the periodic correction:
 [DVectors,IntCurves] = cellfun(@(X,Y) GetPeriodCorrection_DH(X,Y,...
    D,n,Window,SSpacing),SVecs,EVecs,'UniformOutput',0);

%Plot the minima of the axially periodic potentials: 
NumSpirals = 7; %Number of spiral strips
GlobalMinimaPlots(NumSpirals,IntCurves,DVectors,D)
%}
end

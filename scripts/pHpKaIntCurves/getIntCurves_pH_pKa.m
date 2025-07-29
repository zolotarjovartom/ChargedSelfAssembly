% DESCRIPTION: Calculates pairwise interactions between spiral strips as a
% function of pH and residue pKa.

function getIntCurves_pH_pKa
%Specify the collagen type: 
ColType = "II";

%Specify the trimer type: 
TrimerType = "A1|A1|A1";

%Specify the species that you want to perform the calculation for: 
Species = 'BosTaurus';

%Specify the path to raw data files: 
FPathS = '/Volumes/CollagenDrive/Effective Stagger Vectors/';
FPathRes = '/Volumes/CollagenDrive/Bulk Residue Data/';
FPathSFull = strcat(FPathS,'EffDVecs_Mammals_Type',ColType,'|',TrimerType,'.mat');
FPathResFull = strcat(FPathRes,'BulkResData_Type',ColType,'|',TrimerType,'_Mammals.mat'); 

%Add the "methods" path: 
addpath('/Users/artomzolotarjov/Documents/ChargedGitHub/scripts/methods/');

%Specify the path for saving the output file: 
FPathpKa = '/Users/artomzolotarjov/Documents/ChargedGitHub/scripts/constructpHIntCurves/';

%Load the relevant data: 
SVecs = load(FPathSFull).S;
BulkData = load(FPathResFull).s1;

%Create a structure to store output data, generate fieldnames: 
[~,SpeciesNames] = GenStruct(BulkData);

%Delete short sequences:
Lc = 684*4; %Minimal admissable molecular length
BulkData = DeleteShortSeqs(BulkData,Lc);

%Extract the relevant information for further calculations: 
CMat = BulkData.ContactMat; StripCoors = {BulkData.StripAACoors}; 
StripInts = {BulkData.StripAAInts}; 

%Run the calculation: 
pKaGrid(SpeciesNames,Species,SVecs,StripCoors,StripInts,CMat,{ColType,TrimerType},FPathpKa);
end

%SUBROUTINES: 
function [BulkData] = DeleteShortSeqs(BulkData,Lc)
%Extract bulk data:
CaData = {BulkData.CaData};

%For each interaction curve determine molecular length: 
MolL = cellfun(@(X) max(X(:,2))-min(X(:,2)),CaData);

%Impose a limitation on what molecular lengths are admissable:
TF = MolL < Lc;

%Remove short sequence data from BulkData data set: 
BulkData(TF) = [];
Species = BulkData.SeqNamesPaired;

%Amend the list of species:
for fIdx = 1 : numel(BulkData)
    BulkData(fIdx).SeqNamesPaired = Species(~TF);
end
end
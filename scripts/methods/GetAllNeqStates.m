% DESCRIPTION: Construct all near-equilibrium states microstates: 
function [StaggerMat,AzMat] = GetAllNeqStates(v,SSpacing)
%Define the stagger states:
SStates = [1,2,3];
NumAxStates = numel(SStates)^(v-1);

%Construct all possible stagger combinations: 
OutAxStates = cell(1,v-1);
[OutAxStates{1:v-1}] = ndgrid(SStates);
StaggerMat = zeros(NumAxStates,v-1);

for idx = 1:NumAxStates
    StaggerMat(idx,:) = cellfun(@(X) X(idx),OutAxStates);
end

%Construct near equilibrium azimuthal states: 
TStates = GetNearEqbmAngles(v,SSpacing); %Discrete set of near equilibrium azimuthal angles
NumAzimuthStates = numel(TStates)^(v);
OutAzimuthStates = cell(1,v);
[OutAzimuthStates{1:v}] = ndgrid(TStates);

%Create a matrix of near eqbm azimuthal configurations 
AzMat = zeros(NumAzimuthStates,v);

for idx2 = 1:NumAzimuthStates
    AzMat(idx2,:) = cellfun(@(Y) Y(idx2), OutAzimuthStates);
end

end

%SUBROUTINES: 

% Outputs the azimuthal values from which near-equilibrium states will be
% constructred
function [TStates] = GetNearEqbmAngles(v,SSpacing)
EPair = EnergyPairwise_Empirical(1,[1,1;1,1+SSpacing],SSpacing); %Get the azimuthal component of the pairwise interaction energy
EAz = @(az1,az2) EPair(az1,az2);

%Maximise the azimuthal function w.r.t. interior angle of the polygon:
az = linspace(0,2*pi,360*3);
psi = ((v-2)*pi)/v;
[~,locs] = findpeaks(EAz(az,az-psi));

%Output the possible azimuthal angle values for near-eqbm states:
TStates = az(locs);
end

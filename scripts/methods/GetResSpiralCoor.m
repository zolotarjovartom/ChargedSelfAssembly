% DESCRIPTION: outputs spatial coordinates and amino acid names of the
% residues that comprise a selected pair of interacting strips. The spatial
% coordinates of the residues are modified in a way where they form an
% uniterrupted interface in the z-theta plane.

function [StripAACoors] = GetResSpiralCoor(SpiralPairs,SpiralType,StripAACoors)

%Define the z-axis intercepts of the spiral strips:
if SpiralType == 1 %Based on Ca coordinates
    H = 1955.7; %pitch of the spiral strips
    cVals = [-150,-450,-700,-1000,-1250,-1550,-1800];
elseif SpiralType == 2 
    H = 778.4810; %Pitch of the helical strips AMINO    
    cVals= [-40,-110,-190,-270,-350,-430,-510,-590,-670,-750];
end

%Define a convenient parameter:
k = H/(2*pi);

%Adjust the axial coordinates to both start from 0: 
StripAACoors{1}(:,2) = StripAACoors{1}(:,2) - min(StripAACoors{1}(:,2));
StripAACoors{2}(:,2) = StripAACoors{2}(:,2) - min(StripAACoors{2}(:,2));

%Adjust azimuthal coordinates to reflect the translation above: 
StripAACoors{1}(:,1) = mod(StripAACoors{1}(:,1) -  min(StripAACoors{1}(:,2))/k,2*pi);

%Translate the spirals so that they are on top of one another:
SpiralCentreLines = -cVals./k;
ThetaDiff = SpiralCentreLines(SpiralPairs(1)) - SpiralCentreLines(SpiralPairs(2));
StripAACoors{2} = [mod(StripAACoors{2}(:,1) + ThetaDiff - min(StripAACoors{2}(:,2))/k, 2*pi),...
    StripAACoors{2}(:,2:3)];
end
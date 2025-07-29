% DESCRIPTION: Determine which spiral strips contribute to the most stable robust cluster

function [r,c] = GetSpiralStrips(X)
%Determine the number of vertices: 
v = numel(X);

%Next two parameters are the same a sin EnergyPairwise_Empirical
t_p = 0.5; t_d = (2*pi)/7;

%Define the vector containing the internal angles of the polygon: 
psi = ones(1,v).*(-((v - 2)*pi)/v);

%Generate the remaining interacting strips based on the azimuthal constraint
Y = X + psi;

%Map the azimuthal coordinates to integer designations for strips:
StripPairs = zeros(2,v);

%Cyclical cluster:
StripPairs(:,1:v) = [mod(round((X(1:v)-t_p)/t_d),7)+1; mod(round((Y(1:v)-t_p)/t_d),7)+1];
InteractPairs = [1:v;[2:v,1]];

r = StripPairs(1,InteractPairs(1,:));
c = StripPairs(2,InteractPairs(2,:));
end

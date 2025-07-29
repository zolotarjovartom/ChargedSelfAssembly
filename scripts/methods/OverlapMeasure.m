% DESCRIPTION: signed measure of non-overlap between two ionisation intervals.
% The overlap measure is defined as follows: 
% mu = (sign(pKaC(A)-pKaC(B))*(2*delta_pH - |W_{A} n W_{B}|))/2*delta_pH, where n denotes the
% intersection of two ionisation intervals W_{A} and W_{B}. delta_pH
% denotes the half-width of the ionisation intervals (we assume all of them
% are the same width).

% INPUTS: 

% pKaC_A, pKaC_B - the activity adjusted pKa value corresponding to a pH
%with 50% dissociation 

% W_A, W_B - nx2 vectors containing the ionisation intervals for which one
%wants to calculate the signed overlap measure 

% ASSUMPTIONS: 

% The ionisation window along the pH scale is of the same length for all
% residues. 

function mu = OverlapMeasure(pKaC_A,pKaC_B,W_A,W_B)
%Pre-allocate the overlap measure: 
mu = ones(numel(pKaC_A),1); %1 corresponds to complete non-overlap
 
%Identify non-overlapping intervals: 
IsEmpty = or(W_B(:,1) > W_A(:,2),W_A(:,1) > W_B(:,2));

%Identify the ionisation window width: 
width = diff(W_A(1,:)); %Just pick the first element

%Define the overlap measure
if any(~IsEmpty,1)
    nonOverlap = width - (min([W_A(~IsEmpty,2),W_B(~IsEmpty,2)],[],2) - max([W_A(~IsEmpty,1),W_B(~IsEmpty,1)],...
        [],2));
    mu(~IsEmpty) = sign(pKaC_A(~IsEmpty)-pKaC_B(~IsEmpty) + 1.*(abs(pKaC_A-pKaC_B)<1e-4) )...
        .*nonOverlap.*width^(-1);
end
end
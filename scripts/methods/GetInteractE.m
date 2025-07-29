% DESCRIPTION: Given a list of interacting residues, computes the energy of
% interaction at a given stagger

function [InteractionScore] = GetInteractE(AAInts,InteractMat,IntRes)
%Identify the interacting pairs of residues:
StaticIdx = IntRes(:,1); StaggerIdx = IntRes(:,2);
InteractionPairs = sort([AAInts{1}(StaticIdx),AAInts{2}(StaggerIdx)],2,'descend');   

%Compute the interaction energy: 
InteractsToSelect = sub2ind(size(InteractMat),InteractionPairs(:,1),InteractionPairs(:,2));
InteractionScore = sum(InteractMat(InteractsToSelect));
end
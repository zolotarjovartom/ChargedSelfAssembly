% DESCRIPTION: the goal of this script is to redefine pKa values of Arg,
% Asp, Glu and Lys such that the ionisation windows of these residues are
% disjoint intervals.

% INPUTS: 

% IonOrder - order of residue ionisation. The options are: 

% 1 - Glu, Asp, Arg Lys 

% 2 - Asp, Glu, Lys, Arg 

% sf - scale factor, which adds additional separation between ionisation
%intervals.

% alpha_bounds - dissociation fraction which defines the ionisation
%interval. A 2 element vector containing the lower and upper bounds of the
%ionisation interval. 

% NOTE: 

% It is vital that the order of pKa values is preserved. The first 4
% residues needs to be Arg,Asp,Glu,Lys

function [ref_pKa,pKa_shift_out] = Modify_pKa(IonOrder,sf,alpha_bounds,I,...
    stepSize,eps_r,T)
%Function for computing midpoints for intervals:
getMidPt = @(X) X(:,1) + diff(X,1,2)*0.5;

%List of the relevant pKa values:
pKa = [12.3, 3.9, 4.3, 10.4, 9.8, 8.6, 6.5]'; %Arg,Asp,Glu,Lys == 2,4,7,12

%Residue type: 
ResType = {'B','A','A','B','A','A','B'}';

%Vector of charges upon ionisation:
qvec = [1, -1, -1, 1, -1, -1, 1]';

%Compute the width of each dissociation curve: 
alpharange = repmat(alpha_bounds,numel(pKa),1);
[pH_windows,f_act] = get_pHInterval(alpharange,ResType,qvec,eps_r,T,I,pKa);
delta_pH = diff(pH_windows(1,:));

%Identify interval midpoints: 
midPts = getMidPt(pH_windows);

%Update the pKa values s.t. Arg,Asp,Glu,Lys all have pH windows as disjoint
%intervals
pkaMin = min(midPts(end-2:end)); %Minimal pKa from minor ionisable residues
pkaMax = max(midPts(end-2:end)); %Maximal "                               "

%Pick the ordering for ionisation of the first FOUR residues: 
resTypeMod = ResType(1:4); 
IsAcid = matches(resTypeMod,'A')';

%Option 1: ionisation order is Glu, Asp, Arg Lys 
if IonOrder == 1
    ionOrder = [1,-1,-2,2]; %Arg Asp, Glu, Lys, in this order

%Option 2: ionisation order is Asp, Glu, Lys, Arg 
elseif IonOrder == 2
    ionOrder = [2,-2,-1,1];
end

%Generate the new pKa values:
pkaStart = [pkaMax,pkaMin,pkaMin,pkaMax];
activityCorrection = 1.*f_act(1).*IsAcid - 1*f_act(1).*~IsAcid;
ref_pKa = [ionOrder.*sf.*delta_pH + pkaStart+ activityCorrection ...
    ,pKa(end-2:end)'];

%Compute the translation vector for the pKa values: 
stepSize = delta_pH*stepSize;
pKa_shift_out = 0:stepSize:delta_pH;
end

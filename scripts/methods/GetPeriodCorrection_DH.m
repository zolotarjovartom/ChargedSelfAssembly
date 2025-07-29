% DESCRIPTION: when calculating interactions in the DH model we can safely
% assume that for all staggers the interactions are distinct. This was not
% the case when we operated with rectangular potentials.

function [DVec,EMat] = GetPeriodCorrection_DH(DVec,EMat,D,n,Window,SSpacing)
%Determine the first stagger for which the correction is required: 
idxLast = find(diff(EMat(1,:) ),1,'last'); %Last non-zero interaction (positive stagger)
L_Int = mod(-DVec(idxLast),n*D); %Interaction length 
diffD = DVec - L_Int; diffD(diffD <  0) = inf; [~,FirstStagger] = min(diffD);

%Create a list of points at which the correction will be calculated: 
QueryPts = DVec(FirstStagger : idxLast) - n*D;

%Interpolate (linearly) the value of the energy onto the QueryPts:
if size(EMat,1) > 1
    InterpVals = interp1(DVec,EMat',QueryPts)'; %Each COLUMN is a separate entry
elseif size(EMat,1) == 1 
    InterpVals = interp1(DVec,EMat,QueryPts); %Each COLUMN is a separate entry
end

%Add the corrections to the interaction matrix: 
EMat(:,FirstStagger : idxLast) = EMat(:,FirstStagger : idxLast) + InterpVals;

%Truncate the zeros at the end of the energy matrix:
EMat(:,idxLast + 1 : end) = []; DVec(idxLast + 1 : end) = [];

%Truncate unnecessary data: 
diffD = DVec + L_Int; diffD(diffD <=  0) = inf; [~,LastStagger] = min(diffD);
EMat(:,1 : LastStagger - 1) = [];
DVec(1 : LastStagger - 1) = [];

%Identify the region around zero: 
ls = Window(1)*SSpacing; %Smoothing interval around dz is [dz-ls,dz+ls]
idx_lb = find(DVec < -ls - SSpacing,1,"last"); idx_ub = find(DVec > ls,1,"first"); %Note we add an extra point on the left
zero_idx = find(DVec >= 0,1,"first"); %DVec(zero_idx-1 : zero_idx) changes sign
lZero = {EMat(:,idx_lb : zero_idx - 1), DVec(idx_lb : zero_idx - 1)};
rZero = {EMat(:,zero_idx : idx_ub), DVec(zero_idx : idx_ub) + n*D};

%Convert DVec to be modulo n*D: 
DVec = mod(DVec,n*D);
[DVec,SortIdx] = sort(DVec,'ascend');
EMat = EMat(:,SortIdx);

%Append the area around zero (s.t. the periodic b.c. matches): 
EMat = horzcat(lZero{1},EMat,rZero{1}); DVec = horzcat(lZero{2},DVec,rZero{2});
end

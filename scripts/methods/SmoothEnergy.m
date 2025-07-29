% DESCRIPTION: Uses Gaussian-weighted smoothing to compute the smoothed
% interaction energy

function [DVecEqSpaced,SmoothCurve] = SmoothEnergy(EMat,DVec,beta,k,n,D,spacing,FixFlag)
if FixFlag ~= 1
    error('Check how changes in SmoothEnergy.m affect this script!')
end

%Construct a new stagger vector (equally spaced apart from boundary vicinity): 
DVecEqSpaced = horzcat(-(k(1)+1)*spacing : spacing : -spacing, 0 : spacing : spacing*floor((n*D)/spacing), ...
    n*D + spacing : spacing : n*D + k(1)*spacing);
SmoothCurve = zeros(size(EMat,1),size(DVecEqSpaced, 2)- (2*k(1)));

%Smoothing calculation (for each noisy curve):
for idx = 1 : size(EMat,1)
     CurveToSmooth = EMat(idx,:); %Select the noisy curve: 

     CurveToSmooth = interp1(DVec, CurveToSmooth, DVecEqSpaced); %Interpolate the curve onto a uniformly spaced stagger vector
     C = exp(-beta.*CurveToSmooth); %Exponentiate the curve with weight beta
    
     %Select the smoothing algorithm and undo the transformation above: 
     CSmooth = log(smoothdata(C,"gaussian",k))./-beta;
        
     %Only leave entries with staggers in the interval [-spacing,n*D]
     CSmooth(:,1 : k(1)) = []; CSmooth(:,end-k(1)+1 : end) = []; 
     SmoothCurve(idx,:) = CSmooth;    
end

%Truncate the entries outside the interval [0,n*D]:
DVecEqSpaced(1 : k(1)) = []; DVecEqSpaced(end-k(1)+1 : end) = [];
end

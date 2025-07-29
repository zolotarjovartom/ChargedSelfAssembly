% DESCRIPTION: Serves the same purpose as FindFSStates.m but the code is
% vectorised, so this can accept a long vector of inputs.

%NOTES: 

% Dimensionality of input D: this must be a ROW vector of size dim(SState,1)x1

%NOTATION: 

% FS - fully-staggered, understood to have one well-defined axial period 

% CS - complex-staggered, understood to have multiple periodicities. 


%OUTPUTS: 

% TF_FS, TF_CS - a logical vector, which identifies which rows of the SState
% input matrix are classified as FS or CS.

% AxP_FS, AxP_CS - these vectors specify the axial periodicity of FS/CS
% states. The vector entries in AxP_FS correspond to SState(TF_FS,:) and
% similarly for AxP_CS.

function [TF_FS,TF_CS,AxPFS_out,AxPCS_out] = FindFSStates_Vec(SState,D,n,Thresh,v)
%Convert relative staggers into absolute staggers: 
SiMat = mod(cumsum(SState(:,1:v-1),2),n*D);

%Divide through by axial periodicity: 
KVec = 2:5; %Axial period is defined as (n*D)./KVec(i)
AxPMat = repmat(repelem(KVec',size(SState,1),1).*...
    repmat( (n.*D).^(-1),numel(KVec),1 ), 1, size(SState,2) - 1) ;

%AxPMat = repmat(repelem((n.*D).^(-1).*KVec',size(SiMat,1)),1,size(SiMat,2))

SiInt = repmat(SiMat,numel(KVec),1).*AxPMat;

%Check how far off an integer multiple of the axial periodicity we are: 
MixTF = any(abs(SiInt - round(SiInt,0)) > Thresh,2); %Identify periodicities that do not work 

%Identify the set of FS/CFS candidates:
KVals = repelem(KVec',size(SiMat,1)); KCand = KVals(~MixTF);
DVals = repmat(D,numel(KVec),1); DCand = DVals(~MixTF); %The corresponding D-lengthscale values
CandList = horzcat(zeros(size(SiInt(~MixTF,:),1),1),mod(round(SiInt(~MixTF,:),0),KCand));

%Identify the indices of the candidate state in the INPUT SState:
idxCand = mod(find(~MixTF)-1,size(SState,1)) + 1; 

%For the candidate set, calculate the number of unique values:
NUniq = 1 + sum(diff(sort(CandList,2),1,2) ~= 0, 2);

%Construct containers to store the locations/periodicities:
FS_idx = zeros(size(CandList,1),1); %Fully-staggered
CS_idx = zeros(size(CandList,1),1); %Complex-staggered
AxP_FS = zeros(size(CandList,1),1); AxP_CS = zeros(size(CandList,1),1);

for idx = find(NUniq == KCand)'
    [~,~,ix] = unique(CandList(idx,:));
    C = accumarray(ix,1);
  
    %In case all gaps are occupied the same number of times: FS
    if max(abs(diff(C))) == 0
       FS_idx(idx) = idxCand(idx);
       AxP_FS(idx) = n*DCand(idx)*KCand(idx)^(-1);
    %In case there is a single-gap difference:
    elseif max(abs(diff(C))) == 1
       CS_idx(idx) = idxCand(idx);
       AxP_CS(idx) = n*DCand(idx)*KCand(idx)^(-1);
    end
end

%Truncate the zeros: 
CS_idx(CS_idx == 0) = []; FS_idx(FS_idx == 0) = []; 
AxP_CS(AxP_CS == 0) = []; AxP_FS(AxP_FS == 0) = [];

%Construct output files: 
TF_FS = false(size(SState,1),1); TF_CS = false(size(SState,1),1);
TF_FS(FS_idx) = true; TF_CS(CS_idx) = true;
AxPFS_out = zeros(size(SState,1),1); AxPCS_out = zeros(size(SState,1),1);
AxPFS_out(FS_idx) = AxP_FS; AxPCS_out(CS_idx) = AxP_CS;
end
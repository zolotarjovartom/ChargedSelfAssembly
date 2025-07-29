% DESCRIPTION: Computes the equilibrium statistics at a given set of pH and
% pKa values

function [pOut,outStruc] = ComputeGlobalMinimaV2(IntCurve,SVec,DRange,n,v)
%Define parameters for the calculation of periodic interaction energy curves:
Window = [14,14]; %Smoothing is only used for finding minima  
SSpacing = 0.5;
FixFlag = true; beta = 1 ;

%Parameters for state classification:
Thresh = 0.05;
ProbThresh = 0.99;

%Define parameters for the calculation of minima:
NumSpirals = 7;
SpiralPairs = nchoosek(1:NumSpirals,2); 
SpiralPairs = [SpiralPairs;transpose(reshape(repelem(1:NumSpirals,2),[2,NumSpirals]))];

%Initialise the output containers: 
OutCell = cell(numel(DRange),1);

%Preallocate arrays:
pFS = zeros(numel(DRange),1);
pCS = zeros(numel(DRange),1);
Emin = zeros(numel(DRange),1);

for Didx = 1 : numel(DRange)
    %Compute the periodic correction:
     [DVectors,IntCurves] = cellfun(@(X,Y) GetPeriodCorrection_DH(X,Y,...
         DRange(Didx),n,Window,SSpacing),SVec,IntCurve,'UniformOutput',0);

    %Smooth the interactions:    
    [DVectors_S,IntCurves_S] = cellfun(@(X,Y) SmoothEnergy(X,Y,...
        beta,Window,n,DRange(Didx),SSpacing,FixFlag),IntCurves,DVectors,'Un',0);
   
    %Compute the interaction energy minima:  
    SMatrices = ConstructSMatrix_V2(DVectors_S,IntCurves_S,DVectors,IntCurves,DRange(Didx),...
        SpiralPairs,n,NumSpirals,Window.*SSpacing);

    %Use the entire interaction curve, without robustness limitations
    IMat = GetIMat(IntCurves,DVectors,n,DRange(Didx),NumSpirals,SpiralPairs);
    RobIntCurves = {IMat};

    %Determine the most robust FS state:  
    [NeqEnergies,RobStaggers,rob_Az,~,~]...
        = GetEBounds_NoRob(RobIntCurves,{SMatrices},v,SSpacing,n,DRange(Didx));
    
    %Compute the (sorted) cumulative probability:  
    [NeqSort,SortIdx] = sort(NeqEnergies);
    cumProb = cumsum(exp(-1.*NeqSort))./sum(exp(-1.*NeqEnergies));
   
    %Identify the number of (sorted) states required to reach threshold
    %cumulative probability: 
    ThreshIdx = find(cumProb>ProbThresh,1,'first');    
   
    %Classify the NEqS:
    [TF_FS,TF_CS,~,~] = FindFSStates_Vec(RobStaggers(SortIdx(1:ThreshIdx),:),...
        repmat(DRange(Didx),ThreshIdx,1),n,Thresh,v);

    %Compute the finite temperature probabilities of each state: 
    NeqSortTrunc = NeqSort(1:ThreshIdx);
    ZTrunc = sum(exp(-1.*NeqSortTrunc));
    pFS(Didx) = sum(exp(-1.*NeqSortTrunc(TF_FS)))/ZTrunc;
    pCS(Didx) = sum(exp(-1.*NeqSortTrunc(TF_CS)))/ZTrunc;
        
    %Identify the most negative lower bound:
    [LB_min,LB_loc] = min(NeqEnergies);

    %Generate an output structure: 
    OutCell{Didx} = struct('OptStaggers',{RobStaggers(LB_loc,:)},'OptEnergies',{NeqEnergies(LB_loc)},...
            'OptAzStates',{rob_Az(LB_loc,:)});
    Emin(Didx) = LB_min;
end

%Package the probabilities and energies into a structure: 
pOut = struct('pFS',{pFS},'pCS',{pCS},'Emin',{Emin});

%Save the output data:
outStruc = struct('EqStates',{OutCell},'DVals',{DRange});
end
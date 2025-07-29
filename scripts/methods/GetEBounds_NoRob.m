%DESCRIPTION: Compute bounds on the microfibril energy (no robustness)

function [NeqEnergies,rob_Staggers,rob_Az,SStates,Ub_Energies] = GetEBounds_NoRob(RobIntCurves,SMatrices,v,SSpacing,n,D)

%Construct all axial and azimuthal states: 
[StaggerMat,AzMat] = GetAllNeqStates(v,SSpacing);

%Pre-allocate the containers: 
StaggerVals = zeros(7^(v)*3^(v-1),v-1);

%Calculate all the spiral strip compositions and associated terminal interactions:

%Pre-allocate appropriate arrays: 
PairMat = zeros(size(AzMat,1)*v,2);

for azS = 1 : size(AzMat,1)
    %Compute the spiral strip composition: 
    [r,c] = GetSpiralStrips(AzMat(azS,:));
    PairMat( ((azS-1)*v) + 1 : azS*v, : ) = [r;c]'; %Matrix of interacting pairs
end

%Construct the axial state vector: 
Temp = repmat(StaggerMat,1,size(AzMat,1))';
AxVec = [Temp(:)]; %For molecules 1 to v-1

%Convert into linear indices
Idxs = sub2ind([7,7],PairMat(:,1),PairMat(:,2));
Idxs(v:v:end) = []; %Remove indices of terminal interactions
Idxs = repmat(Idxs,size(StaggerMat,1),1);

%Extract the stagger states: 
Staggers = zeros(size(Idxs,1),1);
for minOrd = 1 : numel(SMatrices{1})
    M = AxVec == minOrd; %Identify the interactions of the given order

    %Extract the matrix of stagger states: 
    M_S = SMatrices{1}{minOrd};

    %Populate the stagger matrix: 
    Staggers(M) = M_S(Idxs(M));
end

%Reshape the stagger matrix s.t. each row is a stagger state:
Staggers = reshape(Staggers,size(StaggerVals,2),size(StaggerVals,1))';

%Compute the terminal staggers: 
tStaggers = mod(-sum(Staggers,2),n*D);

%For each unique terminal interaction pair, determine if the terminal stagger is robust: 
IsRob = true(size(tStaggers));

%Identify robust stagger states and corresponding azimuthal states:  
rob_Staggers = [Staggers,tStaggers];
rob_Az = repmat(AzMat,size(StaggerMat,1),1);

%Get azimuthal dependance of the energy: 
EPair = GetAzPenalty;

%Construct mapping vectors: 
mX = 1 : numel(rob_Az);
mY = [2 : numel(rob_Az),1 + (size(rob_Az,1)-1)*v]; mY(v:v:end-1) = 1 + ((1:size(rob_Az,1)-1)-1)*v;

%Compute the azimuthal penalties:
psi = ((v-2)*pi)/v;
TempAz = rob_Az'; rob_AzVec_X = TempAz(:); rob_AzVec_Y = rob_AzVec_X - psi;
AzPenalty = EPair(rob_AzVec_X(mX),rob_AzVec_Y(mY));

%Identify the spiral strips that contribute to the robust states:
Idxs = sub2ind([7,7],PairMat(:,1),PairMat(:,2));
Idxs = repmat(Idxs,size(StaggerMat,1),1);
IsRobVec = repelem(IsRob,v);
IdxsRob = Idxs(IsRobVec);

%Extract the energies of the robust states: 
rob_Staggers_vec = reshape(rob_Staggers',[],1);

%Create a container to store the energies: 
rob_Energies = zeros(numel(rob_Staggers),1);

for Cidx = unique(IdxsRob)'
    %Extract the relevant interaction curve:
    ICurve = RobIntCurves{1}{Cidx};

    %Identify the interactions relevant to this curve:
    M = IdxsRob == Cidx;
    
    %Extract the staggers pertaining to this interaction: 
    curr_Staggers = rob_Staggers_vec(M);
    
    %Compute the unique staggers and identify their associated energies:
    [uniqS,~,idxuniq] =  unique(curr_Staggers);  
    
    %We want to identify the closest value in the interaction curve: 
    edges = [-Inf, mean([ICurve(2,2:end); ICurve(2,1:end-1)]), +Inf];
    I = discretize(uniqS', edges);
    curr_E_vals = ICurve(1,I)';
    rob_Energies(M) = curr_E_vals(idxuniq);
end

%Construct near-equilibrium microfibril energies: 
NeqEnergies = sum(reshape(rob_Energies.*AzPenalty,size(rob_Staggers,2),size(rob_Staggers,1)),1)';

%Construct an upper bound on the microfibril energies:
Ub_Energies = sum(reshape(rob_Energies,size(rob_Staggers,2),size(rob_Staggers,1)),1)';

%Extract the integer stagger states: 
SStates = reshape(AxVec(repelem(IsRob,v-1)),size(rob_Az,2)-1,size(rob_Az,1))';
end
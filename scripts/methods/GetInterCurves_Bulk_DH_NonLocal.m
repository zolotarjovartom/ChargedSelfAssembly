% DESCRIPTION: Computes the electrostatic contribution to the strip
% interaction energy, with relaxed assumption on the locality of
% interactions. For each interacting strip pair i-j we also compute the
% contributions from the interactions between nearest neighbouring strips. 

function [Staggers_pH,Energy_pH] = GetInterCurves_Bulk_DH_NonLocal(NumSpirals,SpiralSet,ContactMat,StripAACoors,StripAAInts,...
    ComputDVec,I,a0,pH,pKaVals,eps_r,T)
%Define pairs of interacting strips: 
SpiralPairs = nchoosek(1:NumSpirals,2);
SpiralPairs = [SpiralPairs;transpose(reshape(repelem(1:NumSpirals,2),[2,NumSpirals]))]; %Include pairs of form n-n
SpiralS = num2cell(SpiralPairs(:,1)); SpiralDyn = num2cell(SpiralPairs(:,2));

%For each pair compute the nearest neighbour interactions: 
LeftStrip = @(sNum) mod(sNum,7) + 1; RightStrip = @(sNum) mod(sNum-2,7) + 1;
LeftPairs = num2cell([LeftStrip(SpiralPairs(:,1)), RightStrip(SpiralPairs(:,2))]);
RightPairs = num2cell([RightStrip(SpiralPairs(:,1)), LeftStrip(SpiralPairs(:,2))]);

%NOTE: by convention, the second COLUMN correspond to the staggered strips

%Set parameters for the calculation: 
CritDist = 7.5; %For non-charged interactions

%This is under the asssumption that we do have constant energy intervals!

%Extract only the left end of each constant energy interval: 
DVecsL = cellfun(@(X) X(1:2:end),ComputDVec,'UniformOutput',0);

%Compute the interaction energies as a function of pH and I:
NumSeeds = numel(pH);
Energy_pH = cell(size(SpiralPairs,1),1); Staggers_pH = cell(size(SpiralPairs,1),1);

for nIdx = 1 : size(SpiralPairs,1)
    [Energy_pH{nIdx},Staggers_pH{nIdx}] = GetNonLocalInts(StripAACoors,[SpiralS{nIdx},LeftPairs{nIdx,1},RightPairs{nIdx,1}], ...
        [SpiralDyn{nIdx},LeftPairs{nIdx,2},RightPairs{nIdx,2}],...
        SpiralSet,StripAAInts,ContactMat,DVecsL,NumSeeds,CritDist,ComputDVec,...
        I,a0,SpiralPairs,pH,pKaVals,eps_r,T);  
end

end

%SUBROUTINES: 

%Compute non-local interactions: 
function [E_out,SVec_out] = GetNonLocalInts(StripAACoors,...
    SpiralNumS,SpiralNumDyn,SpiralSet,StripAAInts,ContactMat,...
    DVecsL,NumSeeds,CritDist,ComputDVec,I,a0,SpiralPairs,pH,pKaVals,eps_r,T)
%Set parameters for the calculation of D-H interactions: 
inv_spacing = 2; % 1/inv_spacing is the step size IN ANGSTROMS in the discretisation of the axial stagger
H = 1955.7; %pitch of the spiral strips 

%Identify the indices that need to be modified in the CP matrix:
pair_vec = [2,4,7,12,19,5,9];
cp_loc = nchoosek(pair_vec,2); cp_loc = [cp_loc;transpose(reshape(repelem(pair_vec,2),[2,numel(pair_vec)]))];

%Modify the CP matrix, such that the CHARGED interaction are 0:
for jIdx = 1 : size(cp_loc,1)
    ContactMat(cp_loc(jIdx,2),cp_loc(jIdx,1)) = 0;  
end

%Compute the strip-strip interactions (including nearest neighbours): 

%Define a container to store the data of the nearest neighbour strips: 
NN_Es = cell(3,1);
NN_SVec = cell(3,1);

for intIdx = 1 : 3
    %Modify the spatial coordinates of the residues to form an interface (azimuthal coordinates match):
    StripPairCoors = {StripAACoors{SpiralNumS(intIdx)},StripAACoors{SpiralNumDyn(intIdx)}};
    [StripPairCoors] = GetResSpiralCoor([SpiralNumS(intIdx),SpiralNumDyn(intIdx)],SpiralSet,StripPairCoors);

    %Extract the residue coordinates + amino acid designations:     
    CaTheta = {StripPairCoors{1}(:,1),StripPairCoors{2}(:,1)};
    CaZ = {StripPairCoors{1}(:,2),StripPairCoors{2}(:,2)};
    CaRho = {StripPairCoors{1}(:,3),StripPairCoors{2}(:,3)};
    StripPairAAInts = {StripAAInts{SpiralNumS(intIdx)},StripAAInts{SpiralNumDyn(intIdx)}};
        
    %Identify the index of the interacting strip pairs: 
    [~,idxPair,~] = intersect(SpiralPairs,sort([SpiralNumS(intIdx),SpiralNumDyn(intIdx)],2),'rows');
  
    %Construct the stagger vector onto which we will interpolate the
    %hydrophobic interactions: 
    FineStagger = ceil(min(ComputDVec{idxPair})) : (1/inv_spacing) : floor(max(ComputDVec{idxPair}));

    %Set the value of parameter that indicates whether or not we are
    %calculating interactions of a nearest neighbour pair: 
    if intIdx == 1 
        NN_TF = false;
    else
        NN_TF = true; 
    end
    
    %Pre-process data for calculation of electrostatic interactions: 
    [dist_vec,z_charge,pKa_vec,RowStatic,RowDyn,V_lc,resType] = PreProcessData_DH_NonLocal(StripPairAAInts,...
        pair_vec,{CaRho,CaTheta,CaZ},inv_spacing,I,H,NN_TF,pKaVals,eps_r,T);
    
    %Check the order of the strips and multiply the stagger vector by -1 if needed:
    if SpiralNumS(intIdx) <= SpiralNumDyn(intIdx)
        staggers = DVecsL{idxPair};

    elseif SpiralNumS(intIdx) > SpiralNumDyn(intIdx)
        staggers = -1.*DVecsL{idxPair};
    end
   
    %Compute the set of all interacting residues for each stagger: 
    IntRes = arrayfun(@(X) GetInteractRes_NonLocal(CaTheta,CaZ,CaRho,X,CritDist,SpiralSet,NN_TF),...
         staggers,'UniformOutput',0);
    
    %Compute the non-electrostatic (hydrophobic) interations, currently pH INDEPENDENT: 
    E_hpob = cell2mat(cellfun(@(X) GetInteractE(StripPairAAInts,ContactMat,X),IntRes,'UniformOutput',0));
    E_hpob = repelem(E_hpob,1,2); %Replicate each energy
    
    %Interpolate the hydrophobic interaction energy onto an
    %interval with the same spacing as the electrostatic interaction

    if SpiralNumS(intIdx) <= SpiralNumDyn(intIdx)
        staggers_full = ComputDVec{idxPair};

    elseif SpiralNumS(intIdx) > SpiralNumDyn(intIdx)
        staggers_full = -1.*ComputDVec{idxPair};
    end 

    E_hpob = interp1(staggers_full,E_hpob,FineStagger);
    
for idx = 1 : NumSeeds 
   
    %On the first pass, determine both the net interaction energy AND 
    % the associated stagger vector SVec_Net  
    if idx == 1    
       %Compute the electrostatic energy (this is pH dependent): 
       [E_El,SVec_El] =  GetAllDH_E(dist_vec,z_charge,pKa_vec,RowStatic,RowDyn,V_lc,pH(idx),I,a0,resType);

       %Compute the stagger vector for the net interaction: 
       [SVec_Net,~,ix] = unique(round(horzcat(FineStagger,SVec_El),5)); %THIS IS THE SOURCE OF DIPLICATION ERROR
       
       %Initialise the container to store the energies: 
       E_mat_pH = zeros(NumSeeds,numel(SVec_Net));
        
       %Add up Hydrophobic + Electrostatic energies and store the net energy: 
       E_mat_pH(idx,:) = accumarray(ix,horzcat(E_hpob,E_El)); 
         
    %For subsequent passes, the stagger vector remains the same, so we only
    %need the interaction energy: 
    elseif idx > 1 
        %Compute the electrostatic energy (this is pH dependent): 
        [E_El,~] =  GetAllDH_E(dist_vec,z_charge,pKa_vec,RowStatic,RowDyn,V_lc,pH(idx),I,a0,resType);
        
        %Store the energy: 
        E_mat_pH(idx,:) = accumarray(ix,horzcat(E_hpob,E_El)); 
    end
    
end
%Store the energy matrices:
NN_Es{intIdx} = E_mat_pH;

%Save the stagger vector:
NN_SVec{intIdx} = SVec_Net;
end

%Sum up all the contributing interactions: 
if NumSeeds > 1
    E_out = NN_Es{1} + interp1(NN_SVec{2},NN_Es{2}',NN_SVec{1})' + interp1(NN_SVec{3},NN_Es{3}',NN_SVec{1})';
elseif NumSeeds == 1 
    E_out = NN_Es{1} + interp1(NN_SVec{2},NN_Es{2},NN_SVec{1}) + interp1(NN_SVec{3},NN_Es{3},NN_SVec{1});
end
SVec_out = NN_SVec{1};

end

% DESCRIPTION: Identifies which residues are interacting at a given stagger

% TAKE NOTE: The inputs into pdist2 go in a SPECIFIC order, namely
% DynamicData FIRST and StaticData SECOND. Also note the change in the order
% of the outputs of find(D<CritDist)

function [IntRes] = GetInteractRes_NonLocal(Theta,Z,Rho,stagger,CritDist,SpiralSet,NN_TF)
%Set the parameters for the calculation: 
   if SpiralSet == 1
      H = 1955.7;
   elseif SpiralSet == 2
      H = 778.4810; 
   end

   %Extract the static and dynamic data: 
   StaticData = [Rho{1},Theta{1},Z{1}];
   DynamicData = [Rho{2},Theta{2},Z{2}];
   
   %Parametrise the metric with pitch and stagger: 
   f_cyl = @(Xi,Yi) cylinder_metric(Xi,Yi,H,stagger,NN_TF);
   
   %Compute the distance matrix: 
   D = pdist2(DynamicData,StaticData,f_cyl);
   
   %NOTE: for pdist(X,Y) D(i,j) corresponds to the distance between
   %entry i in X and entry j in Y.

   %Identify the indices of interacting AAs:
   [StaggerIdx,StaticIdx] = find(D<CritDist);
   IntRes = [StaticIdx,StaggerIdx]; 
end

%Define the metric for the calculation: 

%Xi - a single observation, 1xn vector 

%Yi - multiple observations, m1xn matrix

%H - pitch of the helical strips 

%dZ - stagger between the strips 

%NN_TF - chooses the appropriate metric, depending on whether we are
%calculating nearest neighbour interactions or not.

%1st column is the radius, 2nd column is the azimuth and the 3rd column is
%the axial coordinate

%In order to use the metric with pdist2, pass the parameters H and dZ prior
%to passing the function handle to pdist2. 

function [d_cyl] = cylinder_metric(Xi,Yi,H,dZ,NN_TF)
   if NN_TF == false 
     d_cyl = sqrt(Xi(1).^(2) + Yi(:,1).^(2) ...
        - 2.*Xi(1).*Yi(:,1).*cos(Xi(2) - Yi(:,2) - (dZ.*2.*pi)./H) + (Xi(3)-Yi(:,3)-dZ).^(2));
    
   elseif NN_TF == true 
   %Define the function to compute the radial distance between residues: 
    alpha = (2*pi)/7;
    dp = (Xi(1) + Yi(:,1)) - 2*cos(alpha).*min([repmat(Xi(1),size(Yi,1),1),Yi(:,1)],[],2);
    dn = sqrt(abs(Xi(1)-Yi(:,1)).^(2) + dp.^(2) - 2.*abs(Xi(1)-Yi(:,1)).*dp.*cos(alpha));
    
    %Define the Euclidean metric, corrected for the NN interactions: 
    d_cyl = sqrt(Xi(1).^(2) + (Xi(1)+dn).^(2) ...
        - 2.*Xi(1).*(Xi(1)+dn).*cos(Xi(2) - Yi(:,2) - (dZ.*2.*pi)./H) + (Xi(3)-Yi(:,3)-dZ).^(2));
   end

end

%Pre-process data for non-local DH calculations: 
function [dist_vec,z_charge,pKa,RowStatic,RowDyn,V_lc,resType] = ...
    PreProcessData_DH_NonLocal(AAInts,ResIdx,CylCors,inv_spacing,I,H,NN_TF,pKaVals,...
    eps_r,T)
%Convert the pKa values into a cell of strings: 
pKaVals = strtrim(cellstr(num2str(pKaVals')));

%Find the coordinates of the charged residues: 
TF_Charge = cellfun(@(X) ismember(X,ResIdx),AAInts,'UniformOutput',0);
CylCors_pm = cellfun(@(Rho,Th,Z,Idx) [Rho(Idx), Th(Idx), Z(Idx)],...
    CylCors{1},CylCors{2},CylCors{3},TF_Charge,'UniformOutput',0);

%Around each stagger value, construct a stagger vector: 
V_lc = GenStaggerVecs(inv_spacing,CylCors_pm{1}(:,3) - CylCors_pm{2}(:,3)',I,eps_r,T);

%Associate each stagger vector with appropriate residue charge and pKa: 
pKa = cell(2,1); z_charge = cell(2,1); resType = cell(2,1);
for cIdx = 1 : numel(AAInts)
    %Convert into integer representation: 
    String_Rep = strtrim( cellstr( int2str( AAInts{cIdx}(TF_Charge{cIdx}) ) ) );

    %Correspond with pKa and charge values: 
    z_charge{cIdx} = str2double(replace(String_Rep, {'2';'4';'7';'12';'19';'5';'9'}, {'1';'-1';'-1';'1';'-1';'-1';'1'})); 
    pKa{cIdx} = str2double(replace(String_Rep, {'2';'4';'7';'12';'19';'5';'9'}, pKaVals)); %10.53 Lys; %12.48 Arg

    %Identify whether the residue is basic or acidic:
    resType{cIdx} = replace(String_Rep,{'2';'4';'7';'12';'19';'5';'9'},{'B';'A';'A';'B';'A';'A';'B'});
end 
if NN_TF == false 
    %Define the Euclidean metric in cylindrical coordinates:
    d_cyl = @(ri,rj,ti,tj,zi,zj,dZ,H) sqrt(ri.^(2) + rj.^(2) -2.*ri.*rj.*cos(ti - tj - (dZ.*2.*pi)./H) + (zi-zj-dZ).^(2));
elseif NN_TF == true 
    %Define the function to compute the radial distance between residues: 
    alpha = (2*pi)/7;
    dp = @(ri,rj) (ri + rj) - 2*cos(alpha).*min([ri,rj],[],2);
    dn = @(ri,rj) sqrt(abs(ri-rj).^(2) + dp(ri,rj).^(2) - 2.*abs(ri-rj).*dp(ri,rj).*cos(alpha));
    
    %Define the Euclidean metric, corrected for the NN interactions: 
    d_cyl = @(ri,rj,ti,tj,zi,zj,dZ,H) sqrt(ri.^(2) + (ri+dn(ri,rj)).^(2) ...
        - 2.*ri.*(ri+dn(ri,rj)).*cos(ti - tj - (dZ.*2.*pi)./H) + (zi-zj-dZ).^(2));
end

%Obtain the row indices for CylCors{1} and CylCors{2} that correspond to
%each row of V_lc
NumStatic = size(CylCors_pm{1},1);
NumDyn = size(CylCors_pm{2},1);
NumStaggers = size(V_lc,2); %Number of elements in each stagger vector in V_lc
[RowStatic,RowDyn] = ind2sub([NumStatic,NumDyn], 1 : NumStatic*NumDyn);

%Compute all residue distances: 
%Recall that for each pair of interacting residues we have NumStaggers staggers
RowStatic = repelem(RowStatic,NumStaggers)';
RowDyn = repelem(RowDyn,NumStaggers)';
dist_vec = d_cyl(CylCors_pm{1}(RowStatic,1),CylCors_pm{2}(RowDyn,1),CylCors_pm{1}(RowStatic,2),CylCors_pm{2}(RowDyn,2),...
    CylCors_pm{1}(RowStatic,3),CylCors_pm{2}(RowDyn,3),reshape(V_lc.',[],1),H);
end

%SUBROUTINES: 

% Round to the nearest lattice point:

% DESCRIPTION: suppose we choose to discretise the axial stagger on the
% interval [-L,L], where L is an integer. Suppose this discretisation has 
% step size 1/n, for some integer n. This script rounds all values of Z to
% the nearest discretisation point. 

% inv_spacing - see above

% ZVals - an nxm matrix, defined as [ZVals]ij = Zi-Zj, where Zi and Zj are the
% axial coordinates of residues on strip i and j respectively. 

% I - ionic strength in mM


%OUTPUTS: 

% V_lc - a matrix of size(ZVals,1)xsize(ZVals,2) ROWS. We generate a 
% linearly spaced vector of stagger values ±2 Debye lengths away from a
% each entry in the matrix ZVals. The row index of V_lc is the linear index
% of ZVals. 

function V_lc = GenStaggerVecs(inv_spacing,ZVals,I,eps_r,T)
%Round the staggers to the nearest lattice point:
ZVals = ceil(ZVals.*inv_spacing)./inv_spacing;

%Compute the Debye Length: 

%Define the relevant constants: 
e_charge = 1.602177e-19; %elementary charge
eps_0 = 8.854188e-12; %vaccuum permittivity
k_b = 1.380649e-23; %Boltzmann constant
NA = 6.022141e23; %Avogadro's number

%Compute the Bjerrum length:
l_bjer = e_charge^(2)/(4*pi*eps_0*eps_r*k_b*T); %Bjerrum length in meters

%Debye length:
l_d = ceil(1e10*(1./sqrt(8*pi*l_bjer*NA*I))*inv_spacing)/inv_spacing; %In angstroms and rounded to the nearest lattice point

%Generate linearly spaced stagger vectors for each interaction: 
p = 2; %Number of Debye lengths away form each stagger 
NPoints = p*2*inv_spacing*l_d+1;
x = linspace(0,1,NPoints);
V_lc = ZVals(:)-p*l_d + x.*(2*l_d*p);
end
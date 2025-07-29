% DESCRIPTION: Draws a phase diagram for collagen self-assembly as a
% function of average residue effective charge and mu measure (ionisation
% order & degree of co-ionisation)

function CollagenPhaseDiagram(EqbmData)
%Define variables: 
aRange = [0.05,0.95]; %Defines the width of the ionisation window

%Extract the data:
EqbmData = EqbmData.S; 
fields = fieldnames(EqbmData); TF_mu = contains(fields,'muMeasure'); muFields = fields(TF_mu);
pHMat = EqbmData.pHVec; pKaMat = EqbmData.pKaMat; muVec = EqbmData.muVec;

%Loop over every value of mu measure: 
AclassCell = {}; BclassCell = {};
aNetA = cell(numel(muFields),1); aNetB = cell(numel(muFields),1);

for muIdx = 1 : numel(muFields)
    [tempAClass,tempBClass,aNetA{muIdx},aNetB{muIdx}] ...
        = ClassifyCollagenPhases(EqbmData,muIdx,muFields,pHMat,pKaMat);
    AclassCell = [AclassCell,tempAClass];
    BclassCell = [BclassCell,tempBClass];
end

%Modify the data to be compatible with gscatter: 
figure(3)
tiledlayout(1,2);

%Acidic residues:
nexttile;
hold on;
XDat = cell2mat(aNetA')';
numElem = cellfun(@(X) numel(X),aNetA);
YDat = repelem(muVec,numElem);

%Add a boundary around each cluster: 
uniqCats = unique(AclassCell); %Draw boundaries based on the unique categories in 
%the acid residue classification (the picture is cleaner and can be
%modified if needed).
boundCoors = cell(3,1);
for idx = 1 : numel(uniqCats)
    TF = contains(AclassCell,uniqCats{idx});
    boundX = XDat(TF); boundY = YDat(TF);
    idxB = boundary(boundX,boundY,0.85);
    boundCoors{idx} = [boundX(idxB),boundY(idxB)];
end

%For each boundary, modify the coordinates as appropriate: 
bVec = 1 : 3;
distUB = 0.075;
newBound = cell(3,1);
cols = [0.6350 0.0780 0.1840; 0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250]; counter = 1;
for b = [1,3]
    %Current boundary:
    cB = boundCoors{b};
    
    %Remaining boundary points: 
    P = cell2mat({boundCoors{bVec(bVec~=b)}}');
    
    %Modify each boundary point as appropriate:
    for idx = 1 : size(cB,1)
        [k,dist] = dsearchn(P,cB(idx,:));
        if dist < distUB
           cB(idx,:) = (P(k,:) + cB(idx,:))*0.5;
        end       
    end

    %Smooth portion of the boundary coordinates:
    
    TF1 = and(cB(:,2) < 0.95,cB(:,2) > -0.95);
    TF2 = and(cB(:,1) < 0.95,cB(:,1) > 0.1);
    TFNet = and(TF1,TF2);
    cBMod = smoothdata(cB(TFNet,:),'gaussian',[2,2]);
    cB(TFNet,:) = cBMod;
    
    %Plot the boundary, colour-in:

    fill(cB(:,1),cB(:,2),cols(counter,:),'FaceAlpha',0.5);
    newBound{b} = cB;
    counter = counter + 1;
  
end

%Fill out the rest of the boundary: 
commonB = [newBound{1};newBound{3}];
fillout(commonB(:,1),commonB(:,2),[0.2,0.8,-0.9,0.9],cols(3,:))

%Plot the boundaries: 
lw = 1.5;
plot(newBound{1}(:,1),newBound{1}(:,2),'k-','LineWidth',lw)
plot(newBound{3}(:,1),newBound{3}(:,2),'k-','LineWidth',lw)
plot([0.44,0.93],[1,1],'k-','LineWidth',lw);

%Add the scattered points:

g = gscatter(XDat,YDat,AclassCell','kkk',{'|', 'pentagram', 'x',},6);
for n = 1:length(g)
  set(g(n), 'MarkerFaceColor', 'k');
end
xlabel('$z^{\textrm{ave}}(\textrm{Glu},\textrm{Asp})$',...
    'Interpreter','latex','FontSize',15)
ylabel('$\mu\left(W_{\textrm{Glu}},W_{\textrm{Asp}}\right)$',...
    'Interpreter','latex','FontSize',15)
title('(A)',"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex','FontSize',20,'TitleHorizontalAlignment','left');
xlim([0,1]);
ylim([-1,1]);
grid on;

%Add the curve that splits the phase space into single-coionisation: 
num = aRange(2)*(1-aRange(1));
denom = aRange(1)*(1-aRange(2));
dpH = 0.5*(log10(num/denom));
a_thresh = @(mu) 0.5*(aRange(1) + 1./(1 + 10.^(-2*dpH.*abs(mu)).*(1-aRange(1))*aRange(1)^(-1)));
colAdd = [0.8,0.8,0.8]; lw_add = 3;
plot(a_thresh(-1:0.01:1),-1:0.01:1,'k--','LineWidth',lw_add,'Color',colAdd);

%Compute the mu measure in alanine pentapeptides:
muRef = @(pkA,pkB) sign(pkA-pkB).*(max([pkA,pkB]) - min(pkA,pkB)).*(2*dpH)^(-1);
muRef_acid = muRef(4.3,3.9); muRef_base = muRef(12.3,10.4);

%Plot the reference mu measure value: 
disp(muRef_acid)
disp(muRef_base)
yline(muRef_acid,'LineStyle','-.','LineWidth',lw_add,'Color',colAdd);

%Plot the isoelectronic point (here this is a dummy line) 
%Add a line corresponding to the isoelectronic point of collagen: 
resCounts = load('netResCounts|BosTaurus|TypeII|A1|A1|A1.mat').resCountsOrdered;

%Define the parameters for calculation of ionisation fractions: 
I = 150;  eps_r = 76.5; T = 303.15;
%Define residue types and charges: 
ResType = {'B','A','A','B','A','A','B'}'; qvec = [1, -1, -1, 1, -1, -1, 1]';
z_neut = zeros(size(pKaMat,1),1);
for idx = 1 : size(pKaMat,1)
    pHVec = linspace(min(pHMat{idx}),max(pHMat{idx}),500);
    dissFracs = cell2mat(arrayfun(@(X) GetDissFrac(X,pKaMat(idx,:)',I,qvec,eps_r,T,ResType)',pHVec,'UniformOutput',0)');
    Znet = sum(dissFracs.*repmat(qvec',numel(pHVec),1).*repmat(resCounts,numel(pHVec),1),2);
    [~,minLoc] = min(abs(Znet));
    z_neut(idx) = (dissFracs(minLoc,1) + dissFracs(minLoc,4))/2;
end
plot(NaN,NaN,'LineStyle',':','LineWidth',lw_add,'Color',colAdd);

leg = legend(gca,'Interpreter','latex','Location','northoutside','NumColumns',2);
leglbls = leg.String;
leglbls(end-2) = {'Co-ionisation Boundary'};
leglbls(end-1) = {'$\mu_{\textrm{ref}}(W_{\textrm{R}_{(1)}},W_{\textrm{R}_{(2)}})$'};
leglbls(end) = {'Isoelectronic Point'};
leg.String = leglbls;

hold off;

%Basic residues:
nexttile; 
hold on;
XDat = cell2mat(aNetB')';
numElem = cellfun(@(X) numel(X),aNetB);
YDat = repelem(muVec,numElem);

%NOTE: the markers needs to be adjusted and checked before plotting

%NOTE: for Em threshold set at 0, there is an anomalous point on the boundary

%For illustrative purposes, we take this point out. At finite
%temperatures, due to its proximity to the boundary, we do not expect that
%there microfibrils will form an axially ordered phase, instead the phase
%will not aggregate at all.

XDat(322) = []; YDat(322) = []; BclassCell(322) = [];

%Add a boundary around each cluster: 
boundCoors = cell(3,1);

for idx = 1 : numel(uniqCats)
    TF = contains(BclassCell,uniqCats{idx});
    boundX = XDat(TF); boundY = YDat(TF);
    idxB = boundary(boundX,boundY,0.85);
    boundCoors{idx} = [boundX(idxB),boundY(idxB)];
end

%For each boundary, modify the coordinates as appropriate: 
bVec = 1 : 3;
distUB = 0.075;
newBound = cell(3,1); counter = 1;
for b = [1,3]
    %Current boundary:
    cB = boundCoors{b};
    
    %Remaining boundary points: 
    P = cell2mat({boundCoors{bVec(bVec~=b)}}');
    
    %Modify each boundary point as appropriate:
    for idx = 1 : size(cB,1)
        [k,dist] = dsearchn(P,cB(idx,:));
        if dist < distUB
           cB(idx,:) = (P(k,:) + cB(idx,:))*0.5;
        end       
    end

    %Smooth portion of the boundary coordinates:
    
    TF1 = and(cB(:,2) < 0.95,cB(:,2) > -0.95);
    TF2 = and(cB(:,1) < 0.95,cB(:,1) > 0.1);
    TFNet = and(TF1,TF2);
    cBMod = smoothdata(cB(TFNet,:),'gaussian',[2,2]);
    cB(TFNet,:) = cBMod;
    
    %Edit the boundary points if needed (the procedure above is not perfect):
    if b == 1
        editPts = [0.60187,-0.400548];
        distances = sqrt(sum(bsxfun(@minus, cB, editPts).^2,2));
        %find the smallest distance and use that as an index into B:
        closest = distances==min(distances);
        cB(closest,:) = [0.58,-0.400548];
    end

    %Plot the boundary, colour-in:
    fill(cB(:,1),cB(:,2),cols(counter,:),'FaceAlpha',0.5);
    newBound{b} = cB;
    counter = counter + 1;
  
end

%Fill out the rest of the boundary: 
commonB = [newBound{1};newBound{3}];
fillout(commonB(:,1),commonB(:,2),[0.2,0.8,-0.9,0.9],cols(3,:))

%Plot the boundaries: 
lw = 1.5;
plot(newBound{1}(:,1),newBound{1}(:,2),'k-','LineWidth',lw)
plot(newBound{3}(:,1),newBound{3}(:,2),'k-','LineWidth',lw)
plot([0.43,0.72],[-1,-1],'k-','LineWidth',lw);

%Plot the actual points:
g = gscatter(XDat,YDat,BclassCell','kkk',{'pentagram', 'x','|'},6,'off');

for n = 1:length(g)
  g(n)
  set(g(n), 'MarkerFaceColor', 'k');
end
xlabel('$z^{\textrm{ave}}(\textrm{Arg},\textrm{Lys})$',...
    'Interpreter','latex','FontSize',15)
ylabel('$\mu\left(W_{\textrm{Arg}},W_{\textrm{Lys}}\right)$',...
    'Interpreter','latex','FontSize',15)
title('(B)',"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex','FontSize',20,'xdir','reverse','TitleHorizontalAlignment','left');
xlim([0,1]); grid on;

%Add the curve that splits the phase space into single-coionisation: 
a_thresh = @(mu) 0.5*(aRange(2) + 1./(1 + 10.^(2*dpH.*abs(mu)).*(aRange(2)^(-1)-1)));
plot(a_thresh(-1:0.01:1),-1:0.01:1,'k--','LineWidth',lw_add,'Color',colAdd);

%Add a line corresponding to the isoelectronic point of collagen: 
plot(z_neut,muVec,'LineStyle',':','LineWidth',lw_add,'Color',colAdd);

%Plot the reference mu measure value: 
yline(muRef_base,'LineStyle','-.','LineWidth',lw_add,'Color',colAdd);

hold off;

end

%SUBROUTINES: 

%Classifies the type of aggregate that is predicted at given ionic
%conditions.
function [AclassVec,BclassVec,sum_ad_acid,sum_ad_base] ...
    = ClassifyCollagenPhases(EqbmData,muIdx,muFields,pHMat,pKaMat)
%Define the cluster size: 
N = 5;
eps_tol = 0.05;

%Define the parameters for calculation of ionisation fractions: 
I = 150; 
eps_r = 76.5; %relative permittivity, take value for water at 303.15K
T = 303.15; %Temperature in K 

%Define residue types and charges: 
ResType = {'B','A','A','B','A','A','B'}';
qvec = [1, -1, -1, 1, -1, -1, 1]';

%Extract the data for specific value of mu measure: 
EqbmData = EqbmData.(muFields{muIdx});

%Transform the D variable into gap size:
DRange = EqbmData.DRange; L = 3000; gVar = 5.*DRange - L; %Construct the gap variable

%Extract the pH and pKa vector: 
pHVec = pHMat{muIdx}; pKaVec = pKaMat(muIdx,:);

%Concatenate the probability curves: 
allPCurves = cell2mat(cellfun(@(X) transpose(X.pFS), EqbmData.Probs,'UniformOutput',0));

%Concatenate the energy curves: 
allECurves = cell2mat(cellfun(@(X) transpose(X.Emin), EqbmData.Probs,'UniformOutput',0));

%Compute the energy minima for every pH value: 
[eMins,minLocs] = min(allECurves,[],2);

%Compute the corresponding probabilities: 
rowIdx = transpose(1 : size(allPCurves,1));
pGmin = arrayfun(@(X,Y) allPCurves(X,Y),rowIdx,minLocs);

%Classify the states: 
classVec = cell(1,numel(eMins));

%CONVENTION: 1's - single molecule, 2's - disordered aggregates, 3's - axially-periodic

%Identify the single molecule states: 
eBound = 0; 

%Classification based on infinite microfibrils (just the lowest energy
%state):

TF_SM = eMins./N >= eBound; %Single molecular phase 

%Disordered states: 
TF_DA = false(size(eMins));

%Axially-periodic states: 
TF_AP = false(size(eMins));

%For the states where it is unclear what the state is, identify the lowest
%energy state:
for uncidx = 1 : numel(eMins)
    %Classify based on the staggers of the lowest energy state: 
    OptStaggers = EqbmData.EqStates{uncidx}.EqStates{minLocs(uncidx)}.OptStaggers;
    mathcalD = EqbmData.EqStates{uncidx}.DVals(minLocs(uncidx));
    TF_FS = FindFSStates_Vec(OptStaggers,mathcalD,N,eps_tol,N);

    %If axially periodic: 
    if and(TF_FS,eMins(uncidx) < eBound)
        TF_AP(uncidx) = TF_FS;
    %If disordered:
    elseif and(~TF_FS,eMins(uncidx) < eBound)
        TF_DA(uncidx) = ~TF_FS;
    end
end


%Check that we have accounted for all the data points: 
if ~all(or(or(TF_SM,TF_DA),TF_AP),1)
    warning('Not all states accounted for!')
end

%For the remaning states, 
%Fill in the classification vector: 
classVec(TF_SM) = {'Single Molecule'}; 
classVec(TF_DA) = {'Disordered Aggregate'};
classVec(TF_AP) = {'Axially Periodic'};

%Split the classification vector into acidic and basic components: 
tPoint = find(diff(pHVec) > 1);
AclassVec = classVec(1 : tPoint);
BclassVec = classVec(tPoint + 1 : end);

%Construct the corresponding net disassociation fractions: 
ad = cell2mat(arrayfun(@(X) transpose(GetDissFrac(X,pKaVec',I,qvec,eps_r,T,ResType)),pHVec,'UniformOutput',0)');

%Split up the data: 
ad_acid = ad(1:tPoint,2:3); %Asp, Glu
ad_base = [ad(tPoint+1:end,1),ad(tPoint+1:end,4)]; %Arg, Lys

%Compute the net average disassociation fraction: 
sum_ad_acid = sum(ad_acid,2)'./2;
sum_ad_base = sum(ad_base,2)'./2;

%Extrapolate data to make the boundaries appear better (visually):
sum_ad_acid = [0,sum_ad_acid,1];
AclassVec = ['Single Molecule',AclassVec,'Axially Periodic'];
sum_ad_base = [1,sum_ad_base,0];
BclassVec = ['Axially Periodic',BclassVec,'Single Molecule'];
end

%Compute the disassociation fraction: 
function a_d = GetDissFrac(pH,pKa,I,z,eps_r,T,resType) 
%Set constant(s):
A = (eps_r*T)^(-1.5)*1.82e6;

%Compute the activity coefficient:
I = I*1e-3; %Convert the units of I into M from mM 
f_act = -A.*z.^(2).*(sqrt(I)/(1+sqrt(I)) - 0.2*I);

%Allocate a container for fraction of charged residues: 
a_d = zeros(size(f_act));

%Identify acidic residues: 
TF_acid = matches(resType, 'A');

%Compute the fraction of charged ACIDIC residues: 
a_d(TF_acid) = 1./(1 + 10.^(pKa(TF_acid) - pH + f_act(TF_acid)));

%Compute the fraction of charged BASIC residues: 
a_d(~TF_acid) = 1 - 1./(1 + 10.^(pKa(~TF_acid) - pH - f_act(~TF_acid)));
end 


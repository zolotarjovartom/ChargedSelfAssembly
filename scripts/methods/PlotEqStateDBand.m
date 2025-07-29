% DESCRIPTION: Plots the axial periodicity of the lowest energy state. 

% NOTE: the script is a copy of PaxConservedFeaturesV3.m, but what is being
% plotted on the colorbar is changed. I have not reflected the changes in
% the variable names. 

function PlotEqStateDBand(EqbmData)
%Extract the data:
EqbmData = EqbmData.S;
fields = fieldnames(EqbmData); TF_mu = contains(fields,'muMeasure'); muFields = fields(TF_mu);
pHMat = EqbmData.pHVec; pKaMat = EqbmData.pKaMat; muVec = EqbmData.muVec;

%Loop over different mu values: 
maxP_acid = cell(numel(muFields),1); maxP_base = cell(numel(muFields),1);
a_acid = cell(numel(muFields),1); a_base = cell(numel(muFields),1);
%hold on; 
for muIdx = 1 : numel(muFields)
    [maxP_acid{muIdx},maxP_base{muIdx},a_acid{muIdx},a_base{muIdx}] =...
        plotDbandData(EqbmData,muIdx,muFields,pHMat,pKaMat);
end


%Plot the maximum probability results:

%Acidic residues: 
t = tiledlayout(1,2);

%Containers:
pMat = []; X = []; Y = [];
for idx = 1 : numel(muFields) + 1
    %Plot the corresponding patches:
    if idx < numel(muFields) + 1
        XVec = [0;a_acid{idx}; 1];
        YVec = [-1;muVec;1]; Width = 0.025;
        ub_jidx = numel(maxP_acid{idx}) + 1;
        idxProb = idx;
    else
        XVec = [0;a_acid{end}; 1];
        YVec = [-1;muVec;1]; Width = 0.025;
        ub_jidx = numel(maxP_acid{end}) + 1;
        idxProb = idx - 1;
    end

    for jidx = 1 : ub_jidx
        xRect = [XVec(jidx),XVec(jidx+1),XVec(jidx+1),XVec(jidx)];
        yRect = [YVec(idx)-Width,YVec(idx)-Width,...
            YVec(idx)+Width,YVec(idx)+Width];
    
        X = [X,xRect']; Y = [Y,yRect'];
        %Map the maximum probability onto index vector:
        if jidx < ub_jidx
           pMax = maxP_acid{idxProb}(jidx);
        else 
           pMax = maxP_acid{idxProb}(end);
        end
        pMat = [pMat,pMax];
    end
    
end
h(1) = nexttile(t);
hold on;
%Plot a single patch object:
patch(h(1),X,Y,pMat,'FaceAlpha',1,'EdgeColor','none')

%Overlay a grid over the plot: 
xg = 0:0.1:1;
yg = -1-Width:2*Width:1+Width;
s = [numel(yg) numel(xg)]; % [y x]
hm = mesh(h(1),xg,yg,zeros(s));
hm.FaceColor = 'none'; hm.EdgeColor = [0.5,0.5,0.5]; hm.EdgeAlpha = 0.5;

%Edit the axes: 
FS = 20;
ylabel(h(1),'$\mu(W_{\textrm{Glu}},W_{\textrm{Asp}})$',"Interpreter","latex");
xlabel(h(1),'$z^{\textrm{ave}}(\textrm{Glu},\textrm{Asp})$',"Interpreter","latex");
set(h(1),'TicklabelInterpreter','latex','FontSize',20,'TitleHorizontalAlignment','left');
xlim(h(1),[0,1]); ylim(h(1),[-1,1]); grid on; grid minor; box on;
title("(A)",'Interpreter','latex','FontSize',15)
hold off; 

%Perform the same calculation for basic residues: 
%Containers:

pMat = []; X = []; Y = [];
for idx = 1 : numel(muFields) + 1
    %Plot the corresponding patches:
    if idx < numel(muFields) + 1
       
        XVec = [1;a_base{idx}; 0];
        YVec = [-1;muVec;1]; Width = 0.025;
        ub_jidx = numel(maxP_base{idx}) + 1;
        idxProb = idx;
    else
        XVec = [1;a_base{end}; 0];
        YVec = [-1;muVec;1]; Width = 0.025;
        ub_jidx = numel(maxP_base{end}) + 1;
        idxProb = idx - 1;
    end

    for jidx = 1 : ub_jidx
        xRect = [XVec(jidx),XVec(jidx+1),XVec(jidx+1),XVec(jidx)];
        yRect = [YVec(idx)-Width,YVec(idx)-Width,...
            YVec(idx)+Width,YVec(idx)+Width];
    
        X = [X,xRect']; Y = [Y,yRect'];
        %Map the maximum probability onto index vector:
        if jidx < ub_jidx
           pMax = maxP_base{idxProb}(jidx);
        else 
           pMax = maxP_base{idxProb}(end);
        end
        pMat = [pMat,pMax];
    end
    
end
h(2) = nexttile(t);
hold on;

%Plot a single patch object:
patch(h(2),X,Y,pMat,'FaceAlpha',1,'EdgeColor','none')

%Overlay a grid over the plot: 
xg = 0:0.1:1;
yg = -1-Width:2*Width:1+Width;
s = [numel(yg) numel(xg)]; % [y x]
hm = mesh(h(2),xg,yg,zeros(s));
hm.FaceColor = 'none'; hm.EdgeColor = [0.5,0.5,0.5]; hm.EdgeAlpha = 0.5;

%Edit the axes: 
FS = 20;
ylabel(h(2),'$\mu(W_{\textrm{Arg}},W_{\textrm{Lys}})$',"Interpreter","latex");
xlabel(h(2),'$z^{\textrm{ave}}(\textrm{Arg},\textrm{Lys})$',"Interpreter","latex");
set(h(2),'TicklabelInterpreter','latex','FontSize',20,'XDir','reverse','TitleHorizontalAlignment','left');
xlim(h(2),[0,1]); ylim(h(2),[-1,1]); grid on; grid minor; box on;
title("(B)",'Interpreter','latex','FontSize',15)
hold off; 

% Set colormap and color limits for all subplots
set(h, 'Colormap', parula, 'CLim', [664 669])

% assign color bar to one tile 
cbh = colorbar(h(end),'TickLabelInterpreter','latex'); 

% To position the colorbar as a global colorbar representing
% all tiles, 
cbh.Layout.Tile = 'east'; 

%Add a colorbar:
ylabel(cbh,'$\mathcal{D}\,/\,\textrm{\AA}$','Interpreter','latex',...
    'FontSize',FS)
end

function [AxPFS_acid,AxPFS_base,sum_ad_acid,sum_ad_base] = plotDbandData(EqbmData,muIdx,muFields,pHMat,pKaMat) 
%Define the parameters for calculation of ionisation fractions: 
I = 150; 
eps_r = 76.5; %relative permittivity, take value for water at 303.15K
T = 303.15; %Temperature in K 
aRange = [0.05,0.95]; %Defines the width of the ionisation window
N = 5;
eps_tol = 0.01;
%Define residue types and charges: 
ResType = {'B','A','A','B','A','A','B'}';
qvec = [1, -1, -1, 1, -1, -1, 1]';

%Extract the data for specific value of mu measure: 
EqbmData = EqbmData.(muFields{muIdx});

%Transform the D variable into gap size:
DRange = EqbmData.DRange; L = 3000; gVar = 5.*DRange - L; %Construct the gap variable

%Extract the pH and pKa vector: 
pHVec = pHMat{muIdx}; pKaVec = pKaMat(muIdx,:);

%Concatenate the energy curves: 
allECurves = cell2mat(cellfun(@(X) transpose(X.Emin), EqbmData.Probs,'UniformOutput',0));

%Compute the energy minima for every pH value: 
[eMins,minLocs] = min(allECurves,[],2);

%Identify the single molecular states: 
eBound = 0;
TF_SM = eMins./N >= eBound; %Single molecular phase 

%Disordered states: 
TF_DA = false(size(eMins));

%Construct a vector to store the periodicities: 
AxPFS = zeros(numel(eMins),1);

for uncidx = 1 : numel(eMins)
    %Classify based on the staggers of the lowest energy state: 
    OptStaggers = EqbmData.EqStates{uncidx}.EqStates{minLocs(uncidx)}.OptStaggers;
    mathcalD = EqbmData.EqStates{uncidx}.DVals(minLocs(uncidx));
    [TF_FS,~,AxPFS(uncidx),~] = FindFSStates_Vec(OptStaggers,mathcalD,N,eps_tol,N);
    
     %If axially periodic: 
    if and(~TF_FS,eMins(uncidx) < eBound)
        TF_DA(uncidx) = ~TF_FS;
    end
end

%Keep only the axial periods corresponding to non single molecule phase: 
AxPFS(or(TF_SM,TF_DA)) = NaN;

%Modify the data point that lies on the boundary: 
if muIdx == 9 
    AxPFS(51) = NaN;
end

%Identify the points where we switch over from looking at acidic to basic residues:
tPoint = find(diff(pHVec) > 1);

%Compute the ionisation fractions
a_d = cell2mat(arrayfun(@(X) transpose(GetDissFrac(X,pKaVec',I,qvec,eps_r,T,ResType)),pHVec,'UniformOutput',0)');

%Split up the data: 
a_d_acid = a_d(1:tPoint,2:3); %Asp, Glu
a_d_base = [a_d(tPoint+1:end,1),a_d(tPoint+1:end,4)]; %Arg, Lys

%Compute the average disassociation fraction: 
sum_ad_acid = sum(a_d_acid,2)./2;
sum_ad_base = sum(a_d_base,2)./2;

%Compute the maximum probability for each curve: 
AxPFS_acid = AxPFS(1:tPoint);
AxPFS_base = AxPFS(tPoint+1 : end);
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



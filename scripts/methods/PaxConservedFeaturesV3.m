% DESCRIPTION: in the mu-z^{ave} space, compute the maximum equilibrium
% probability of axially ordered aggregates.

function PaxConservedFeaturesV3(EqbmData)
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
        plotProbCurveData(EqbmData,muIdx,muFields,pHMat,pKaMat,muVec);
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

%Add the curve that splits the phase space into single-coionisation: 
aRange = [0.05,0.95];
num = aRange(2)*(1-aRange(1));
denom = aRange(1)*(1-aRange(2));
dpH = 0.5*(log10(num/denom));
a_thresh = @(mu) 0.5*(aRange(1) + 1./(1 + 10.^(-2*dpH.*abs(mu)).*(1-aRange(1))*aRange(1)^(-1)));
colAdd = [0.8,0.8,0.8]; lw_add = 3;
plot(h(1),a_thresh(-1:0.01:1),-1:0.01:1,'k--','LineWidth',lw_add,'Color',colAdd);

%Overlay a grid over the plot: 
xg = 0:0.1:1;
yg = -1-Width:2*Width:1+Width;
s = [numel(yg) numel(xg)]; % [y x]
hm = mesh(h(1),xg,yg,zeros(s));
hm.FaceColor = 'none'; hm.EdgeColor = [0.5,0.5,0.5]; hm.EdgeAlpha = 0.5;

%Edit the axes: 
ylabel(h(1),'$\mu(W_{\textrm{Glu}},W_{\textrm{Asp}})$',"Interpreter","latex");
xlabel(h(1),'$z^{\textrm{ave}}(\textrm{Glu},\textrm{Asp})$',"Interpreter","latex");
set(h(1),'TicklabelInterpreter','latex','FontSize',20);
xlim(h(1),[0,1]); ylim(h(1),[-1,1]); grid on; grid minor; box on;
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

%Add the curve that splits the phase space into single-coionisation: 
a_thresh = @(mu) 0.5*(aRange(2) + 1./(1 + 10.^(2*dpH.*abs(mu)).*(aRange(2)^(-1)-1)));
cib = plot(h(2),a_thresh(-1:0.01:1),-1:0.01:1,'k--','LineWidth',lw_add,'Color',colAdd);

legend(cib,'Co-ionisation Boundary','Interpreter','latex','FontSize',20,'Location','NorthOutside')
%Edit the axes: 
FS = 20;
ylabel(h(2),'$\mu(W_{\textrm{Arg}},W_{\textrm{Lys}})$',"Interpreter","latex");
xlabel(h(2),'$z^{\textrm{ave}}(\textrm{Arg},\textrm{Lys})$',"Interpreter","latex");
set(h(2),'TicklabelInterpreter','latex','FontSize',20,'XDir','reverse');
xlim(h(2),[0,1]); ylim(h(2),[-1,1]); grid on; grid minor; box on;
hold off; 

% Set colormap and color limits for all subplots
set(h, 'Colormap', hot, 'CLim', [0 1])

% assign color bar to one tile 
cbh = colorbar(h(end),'TickLabelInterpreter','latex'); 

% To position the colorbar as a global colorbar representing
% all tiles, 
cbh.Layout.Tile = 'east'; 

%Add a colorbar:
ylabel(cbh,'$\max\limits_{g}\{\mathcal{P}_{\textrm{period}}(g)\}$','Interpreter','latex',...
    'FontSize',FS)
end

%SUBROUTINES:

function [maxP_acid,maxP_base,sum_ad_acid,sum_ad_base] = plotProbCurveData(EqbmData,muIdx,muFields,pHMat,pKaMat,muVec) 
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
maxP_acid = max(allPCurves(1:tPoint,:),[],2);
maxP_base = max(allPCurves(tPoint + 1:end,:),[],2);
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

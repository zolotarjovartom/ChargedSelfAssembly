% DESCRIPTION: This script plots the values of E_m^{min} and
% P_{ax}(g_{min}) as a function of z^{ave}. 

function PlotEqProb_EqEnergy_V2(EqbmData,muSlice)
%Extract the data:
EqbmData = EqbmData.S; 
fields = fieldnames(EqbmData); TF_mu = contains(fields,'muMeasure'); muFields = fields(TF_mu);
pHMat = EqbmData.pHVec; pKaMat = EqbmData.pKaMat; muVec = EqbmData.muVec;

%Select the mu measure slices: 
N = 5; %Size of the microfibril

%Display the mu measure slices selected:
disp(muFields{muSlice})

%Define the containers: 
pGmin  = cell(numel(muSlice),1); eMins = cell(numel(muSlice),1);
stateClass_acid = cell(numel(muSlice),1); stateClass_base = cell(numel(muSlice),1);
sum_ad_acid = cell(numel(muSlice),1); sum_ad_base = cell(numel(muSlice),1);

%Extract the necessary data for each slice: 
for idx = 1 : numel(muSlice)
    [pGmin{idx},eMins{idx},stateClass_acid{idx},stateClass_base{idx},...
        sum_ad_acid{idx},sum_ad_base{idx}] = ComputeEqbmClassifiers(EqbmData,muSlice(idx),muFields,pHMat,pKaMat);
end


%Define the state colours: 
curveCols = [0.4940 0.1840 0.5560 0.5; ...
            0.9290 0.6940 0.1250 0.5 ;...
            0.6350 0.0780 0.1840 0.5; ...
            [0, 0.4470, 0.7410], 0.5]; 

%Plot the data: 
figure(muSlice)
%This is the global plot Layout:
T = tiledlayout(2,numel(pGmin),"TileSpacing","tight",'TileIndexing','rowmajor');

%The probability curves: 
TFProb = true;
numel(pGmin)
for idx = 1 : numel(pGmin)
    %Create the nested axes inside hte global Layout:
    t = tiledlayout(T,1,2,'TileSpacing','compact');
    bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
    bgAx.Layout.TileSpan = [1 2];
    t.Layout.Tile = idx;
    
    %Define the inputs (acid):
    YIn = pGmin{idx}(1:numel(sum_ad_acid{idx}));
    XIn = sum_ad_acid{idx};
    Classifier = stateClass_acid{idx};
    
    %First axis:
    ax1 = axes(t);
    tile = 1;
    [~,tCoors] = plotCurveWithTranRegion(XIn,YIn,Classifier,curveCols,ax1,tile,TFProb);
      
    %Add the titles:
    titleStr = strcat('$\mu(W_{\textrm{R}_{(1)}},W_{\textrm{R}_{(2)}})=',num2str(muVec(muSlice(idx))),'$');
    title(titleStr,'Interpreter','latex','Position',[1.1,1.05])
    
    %Edit the first axis: 
    FS = 20;
    if idx == 1
        ylabel(ax1,'$\lim\limits_{L_{\textrm{M}}\to\infty}\mathcal{P}_{\textrm{period}}(g_{\textrm{min}})$',...
           'Interpreter','latex','FontSize',FS)
    end
    grid(ax1,'on'); grid(ax1,'minor')
    xlim(ax1,[0,1.0]); ylim(ax1,[-0.01,1.0])
    
    %Define the inputs (base):
    YIn = pGmin{idx}(numel(sum_ad_acid{idx})+1:end);
    XIn = sum_ad_base{idx};
    Classifier = stateClass_base{idx};
    
    %Second axis: 
    ax2 = axes(t);
    tile = 2;
    plotCurveWithTranRegion(XIn,YIn,Classifier,curveCols,ax2,tile,TFProb);
    set(ax2,'XDir','reverse','box','off'); ax2.YAxis.Visible = 'off';  
    grid(ax2,'on');grid(ax2,'minor'); xlim(ax2,[0,1]); ylim(ax2,[-0.01,1.0])

end

%The microfibril energy curves: 
TFProb = false;
for idx = 1 : numel(pGmin)
    %Create the nested axes inside hte global Layout:
    t = tiledlayout(T,1,2,'TileSpacing','compact');
    bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
    bgAx.Layout.TileSpan = [1 2];
    t.Layout.Tile = idx + 1;
    
    %Define the inputs (acid):
    YIn = eMins{idx}(1:numel(sum_ad_acid{idx}))./N;
    XIn = sum_ad_acid{idx};
    Classifier = stateClass_acid{idx};
    
    %First axis:
    ax1 = axes(t);
    tile = 1;
    [tLbls,tCoors] = plotCurveWithTranRegion(XIn,YIn,Classifier,curveCols,ax1,tile,TFProb);
    
    %Edit the first axis: 
    FS = 20;
    xlabel(ax1,'$z^{\textrm{ave}}(\textrm{Glu},\textrm{Asp})$',...
        'Interpreter','latex','FontSize',FS)
    ylabel(ax1,'$\textrm{E}_{\textrm{M}}^{\textrm{min}}\textrm{N}^{-1}\,/\,k_{B}T$',...
        'Interpreter','latex','FontSize',FS)
    grid(ax1,'on'); grid(ax1,'minor')
    xlim(ax1,[0,1]); ylim(ax1,[-15,15])

    %Add a second x-axis: 
    axS = axes(t);   
    axS.XAxisLocation = 'top'; axS.YAxis.Visible = 'off'; axS.Color = 'none'; axS.Box = 'off';
    set(axS,'TickLabelInterpreter','latex','FontSize',20,'box','off');
    xlim(axS,[0,1.0]); xticks(axS,tCoors(:,1)); xticklabels(axS,tLbls');
    xtickangle(axS,0);
    axS.XRuler.TickLabelGapOffset = 1;

    %Define the inputs (base):
    YIn = eMins{idx}(numel(sum_ad_acid{idx})+1:end)./N;
    XIn = sum_ad_base{idx};
    Classifier = stateClass_base{idx};
    
    %Second tile: 
    ax2 = axes(t);
    tile = 2;
    [tLbls,tCoors] = plotCurveWithTranRegion(XIn,YIn,Classifier,curveCols,ax2,tile,TFProb);
    [tCoors,sortIdx] = sortrows(tCoors,1); 
    tLbls = tLbls(sortIdx);
    set(ax2,'XDir','reverse','box','off'); ax2.YAxis.Visible = 'off';  
    grid(ax2,'on');grid(ax2,'minor'); 
    ylim(ax2,[-15,15])
    xlabel(ax2,'$z^{\textrm{ave}}(\textrm{Arg},\textrm{Lys})$',...
        'Interpreter','latex','FontSize',FS);
    
    %Add a second x-axis: 
    axS2 = axes(t); 
    axS2.Layout.Tile = 2;
    axS2.XAxisLocation = 'top'; axS2.YAxis.Visible = 'off'; axS2.Color = 'none'; axS.Box = 'off';
    set(axS2,'TickLabelInterpreter','latex','FontSize',20,'XDir','reverse');
    xlim(axS2,[0,1.0]); xticks(axS2,tCoors(:,1)); xticklabels(axS2,tLbls');
    xtickangle(axS2,0);
    axS2.XRuler.TickLabelGapOffset = 1;
end

end

%SUBROUTINES: 

%Computes P_{ax}(g_{min}) and E_{M}^{min} curves as a function of z^{ave}
function [pGmin,eMins,stateClass_acid,stateClass_base,sum_ad_acid,sum_ad_base] = ...
    ComputeEqbmClassifiers(EqbmData,muIdx,muFields,pHMat,pKaMat)
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

%Compute the probabilities for an infinite microfibril: 

%Identify the single molecule states: 
eBound = 0; 
TF_SM = eMins./N >= eBound; 

%Identify the remaining states and assign probabilities:
pGmin = zeros(size(eMins));

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
        pGmin(uncidx) = 1;
    %If disordered:
    elseif and(~TF_FS,eMins(uncidx) < eBound)
        TF_DA(uncidx) = ~TF_FS;
        pGmin(uncidx) = 0;
    %If single molecule
    elseif TF_SM(uncidx)
        if TF_FS
           pGmin(uncidx) = 1; 
        else
           pGmin(uncidx) = 0;
        end     
    end

end

%Check that we have accounted for all the data points: 
if ~all(or(or(TF_SM,TF_DA),TF_AP),1)
    warning('Not all states accounted for!')
end

%Concatenate the logical arrays: 
stateClass = [TF_SM,TF_DA,TF_AP];
%}

%Construct the net disassociation fractions: 
ad = cell2mat(arrayfun(@(X) transpose(GetDissFrac(X,pKaVec',I,qvec,eps_r,T,ResType)),pHVec,'UniformOutput',0)');

%Split up the disassociation fraction data: 
tPoint = find(diff(pHVec) > 1);
ad_acid = ad(1:tPoint,2:3); %Asp, Glu
ad_base = [ad(tPoint+1:end,1),ad(tPoint+1:end,4)]; %Arg, Lys

%Compute the net average disassociation fraction: 
sum_ad_acid = sum(ad_acid,2)'./2;
sum_ad_base = sum(ad_base,2)'./2;

%Split up the state classification data:
stateClass_acid = stateClass(1:tPoint,:);
stateClass_base = stateClass(tPoint+1:end,:);
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

%Plot the relevant curve indicating the transition region: 
function [tLbls,tCoors] = plotCurveWithTranRegion(XIn,YIn,Classifier,...
    curveCols,ax,tile,TFProb)
%Plot the data:
Mrkrs ={'|','x','pentagram','diamond'}; LW = 8;
hold on; 
ax.Layout.Tile = tile;

%Define a matrix to store transition point coordinate:
tCoors = [];
for classIdx = 1 : size(Classifier,2)
    
    %Extract the section of the curve corresponding to phase classIdx:
    YData = YIn(Classifier(:,classIdx));
    XData = XIn(Classifier(:,classIdx));

    %Check the location of the endpoint:
    stateLoc = find(Classifier(:,classIdx));
    if ~isempty(stateLoc)
     stateLocEnd = stateLoc(end);
        if stateLocEnd ~= size(Classifier(:,classIdx),1)
           %Plot the connecting region between phases and mark the transition point: 
           XMidPt = mean(XIn(stateLocEnd:stateLocEnd+1)); YMidPt = mean(YIn(stateLocEnd:stateLocEnd+1));
           X = [XData(end),XMidPt]; Y = [YData(end),YMidPt];
           p1 = patch([X,flip(X)],[Y,flip(Y)],curveCols(classIdx,1:3));
           if ~isempty(p1)
            p1.EdgeColor=curveCols(classIdx,1:3);
            p1.FaceColor='none';
            p1.EdgeAlpha=curveCols(classIdx,end);
            p1.LineWidth=LW;
            p1.LineJoin='round'; % to make it similar to plot's default
           end

           colTF = Classifier(stateLocEnd+1,:)';
           X = [XMidPt,XIn(stateLocEnd+1)]; Y = [YMidPt,YIn(stateLocEnd+1)];
           p2 = patch([X,flip(X)],[Y,flip(Y)],curveCols(classIdx,1:3));
           if ~isempty(p2)
            p2.EdgeColor=curveCols(colTF,1:3);
            p2.FaceColor='none';
            p2.EdgeAlpha=curveCols(colTF,end);
            p2.LineWidth=LW;
            p2.LineJoin='round'; % to make it similar to plot's default
           end

           %Store the coordinates: 
           tCoors = [tCoors;XMidPt,YMidPt,classIdx,find(colTF)];
        
        end

    end

    %Plot the curve portion: (does not preserve transperency)
    %plot(ax,XData,YData,'Color',curveCols(classIdx,:),...
    %    'LineWidth',LW)

    p = patch(ax,[XData flip(XData)],[YData' flip(YData')],curveCols(classIdx,1:3));
    if ~isempty(p)
        p.EdgeColor=curveCols(classIdx,1:3);
        p.FaceColor='none';
        p.EdgeAlpha=curveCols(classIdx,end);
        p.LineWidth=LW;
        p.LineJoin='round'; % to make it similar to plot's default
    end
    %Add the markers: 
     plot(ax,XData,YData,'Color',curveCols(classIdx,:),...
        'LineStyle','none','LineWidth',2,'Marker',Mrkrs{classIdx},'MarkerSize',6,"MarkerFaceColor",'k',"MarkerEdgeColor",'k')
end
%Add vertical lines:
if ~TFProb
    for lidx = 1 : size(tCoors,1)
        plot(ax,[tCoors(lidx,1),tCoors(lidx,1)],[tCoors(lidx,2),15],':','LineWidth',2,'Color','k')
        plot(ax,tCoors(lidx,1),tCoors(lidx,2),'ko','LineWidth',2.5,'MarkerSize',10)
    end
else 
    for lidx = 1 : size(tCoors,1)
        plot(ax,[tCoors(lidx,1),tCoors(lidx,1)],[tCoors(lidx,2),-0.1],':','LineWidth',2,'Color','k')
        plot(ax,tCoors(lidx,1),tCoors(lidx,2),'ko','LineWidth',2.5,'MarkerSize',10)
    end
end
%Label the transitions: 
phaseLbls = {'\vert';'\times';'\star';'\diamond'};
tLbls = strcat('$\textrm{z}^{\textrm{ave}}_{',phaseLbls(tCoors(:,3)),'\to',phaseLbls(tCoors(:,4)),'}$');

%Edit the axes:
set(ax,'TickLabelInterpreter','latex','FontSize',20);

%Add a line to mark the end: 
xline(1,'LineWidth',1.5,'LineStyle','--')
hold off;  

end
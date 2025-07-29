%DESCRIPTION: A more generalised version of PlotEqProb_EqEnergy_pH. The
%idea is to plot the profiles of P_{ax}(g) and min_{S}{E_{M}(g)} around the
%phase transition points and identify distinct mechanisms involved in the
%phase transition.

%This script is designed to work with data sets that involve the \mu
%measure. 

function PlotEqProb_EqEnergy_V3(EqbmData,inputIdx)
%Extract the data:
EqbmData = EqbmData.S; 
fields = fieldnames(EqbmData); TF_mu = contains(fields,'muMeasure'); muFields = fields(TF_mu);

%Select the mu measure slices: 
muSlice = inputIdx;
disp(muFields{inputIdx})

%Extract the necessary data for a given slice: 
[pData,eData,phaseClassifier,gVar,zArg] =...
    ComputeEqbmClassifiers(EqbmData,muSlice,muFields);

%The plot layout:
figure(inputIdx)
T = tiledlayout(3,2,'TileSpacing','compact','TileIndexing','rowmajor');

%Loop over the differen phase transitions: 
colororder({'k','b'}); FS = 15; FSlbl = 25;
counter = 1; 
for pidx = 1 : numel(phaseClassifier)
   
    %Loop over the profiles over different z^{ave} values:
    for zidx = 1 : size(pData{pidx},1)
        ax = nexttile;
        hold on;
        
        %Add normalised energy to the plot: 

        %Determine the normalisation factor: 
        yyaxis(ax,'right')
        if max(eData{pidx}(zidx,:)) > 0 
            NFac = max(abs(eData{pidx}(zidx,:)));
            plot(ax,gVar,eData{pidx}(zidx,:)./NFac,'b-','LineWidth',2.0)
            ylim(ax,[-1,1])
            set(gca,'TickLabelInterpreter','latex','FontSize',FS)
            
        else
            NFac = abs(min(eData{pidx}(zidx,:)));
            plot(ax,gVar,eData{pidx}(zidx,:)./NFac,'b-','LineWidth',2.0)
            ylim([-1,0])
        end
        %Plot the energy threshold: 
        %yline(ax,-1./NFac,'b:','LineWidth',2); % -1 k_{B}T threshold 
        yline(ax,0,'b:','LineWidth',2); % 0 k_{B}T threshold

        %Plot the probability:
        yyaxis(ax,'left')
        plot(ax,gVar, pData{pidx}(zidx,:),'k-',"LineWidth",1.5)
        
        %Plot the probability threshold: 
        yline(ax,0.95,'k--','LineWidth',2);
        ylim(ax,[0,1]); xlim([min(gVar),max(gVar)]);
        grid on; grid minor;
       
        %Add a title: 
        if zidx == 1
            title(strcat(zArg{pidx},'<',phaseClassifier{pidx}),'Interpreter','latex')
        elseif zidx == 2
             title(strcat(zArg{pidx},'>',phaseClassifier{pidx}),'Interpreter','latex')
        end
        
        set(gca,'TickLabelInterpreter','latex','FontSize',FS)
        hold off;
        counter = counter + 1;
    end
end
% X label
xlabel(T,'$g\,/\,\textrm{\AA}$',"Interpreter","latex",'FontSize',FSlbl)

% Left label
ylabel(T,'$\mathcal{P}_{\textrm{period}}(g)$',"Interpreter","latex",'FontSize',FSlbl)
end

%SUBROUTINES: 

%Computes P_{ax}(g_) and min_{S}{E_{M}(g)} curves around the phase
%transition point.

function [pData,eData,tNames,gVar,zArg] = ComputeEqbmClassifiers(EqbmData,muIdx,muFields)
%Define the cluster size: 
N = 5; 
eps_tol = 0.05;
%Extract the data for specific value of mu measure: 
EqbmData = EqbmData.(muFields{muIdx});

%Transform the D variable into gap size:
DRange = EqbmData.DRange; L = 3000; gVar = 5.*DRange - L; %Construct the gap variable

%Concatenate the probability curves: 
allPCurves = cell2mat(cellfun(@(X) transpose(X.pFS), EqbmData.Probs,'UniformOutput',0));

%Concatenate the energy curves: 
allECurves = cell2mat(cellfun(@(X) transpose(X.Emin), EqbmData.Probs,'UniformOutput',0));

%Compute the energy minima for every pH value: 
[eMins,minLocs] = min(allECurves,[],2);

%Compute the corresponding probabilities: 
rowIdx = transpose(1 : size(allPCurves,1));
pGmin = arrayfun(@(X,Y) allPCurves(X,Y),rowIdx,minLocs);

%Identify the single molecule states: 
eBound = 0; pBound = 0.95;

TF_SM = eMins./N >= eBound; 

%Identify the disordered states: 
TF_DA = and(eMins./N < eBound,pGmin <= 1 - pBound);

%Identify the axially-periodic states: 
TF_AP = and(eMins./N < eBound,pGmin >= pBound);

%Identify states that do not correspond to the above classifications
TF_Unclear = and(eMins./N < eBound, and(pGmin > 1-pBound,pGmin < pBound));

%For the states where it is unclear what the state is, identify the lowest
%energy state:
if sum(TF_Unclear) > 0
    for uncidx = find(TF_Unclear)'
        %Classify based on the staggers of the lowest energy state: 
        OptStaggers = EqbmData.EqStates{uncidx}.EqStates{minLocs(uncidx)}.OptStaggers;
        mathcalD = EqbmData.EqStates{uncidx}.DVals(minLocs(uncidx));
        TF_FS = FindFSStates_Vec(OptStaggers,mathcalD,N,eps_tol,N);
        
        %If axially periodic: 
        if TF_FS
            TF_AP(uncidx) = TF_FS;
        %If disordered:
        else
            TF_DA(uncidx) = ~TF_FS;
        end
    end
end

%Check that we have accounted for all the data points: 
if ~all(or(or(TF_SM,TF_DA),TF_AP),1)
    warning('Not all states accounted for!')
end

%Assign: 1 --> SM, 2 --> DA, 3 --> AP
stateClass = zeros(size(eMins));
stateClass(TF_SM) = 1; stateClass(TF_DA) = 2; stateClass(TF_AP) = 3;
%Identify the phase transitions and categorize them:
tPoints = find(diff(stateClass)~=0);
lD = 1; rD = 2; %Specify the distance from the transition point: 
tDoublets = [stateClass(tPoints-lD),stateClass(tPoints+rD)]; %In a given row
%left to right columns indicate a phase transition with increasing z^{ave}

%Identify acidc/basic transitions: 
TF_acid = diff(tDoublets,1,2)>0;
zArg = cell(size(tDoublets,1),1);
zArg(TF_acid) = {'$z^{\textrm{ave}}(\textrm{Glu},\textrm{Asp})'};
zArg(~TF_acid) = {'$z^{\textrm{ave}}(\textrm{Arg},\textrm{Lys})'};

%Classify the phase transitions: 
Mrkrs = {'\vert';'\times';'\star';'\diamond'};
tNames = strcat('z^{\textrm{ave}}_{',Mrkrs(tDoublets(:,1)),'\to',Mrkrs(tDoublets(:,2)),'}$');

%Extract the relevant probability/energy profiles:
eData = cell(numel(tPoints),1); %NOTE: energy is per molecule
pData = cell(numel(tPoints),1);
for idx = 1 : numel(tPoints)
    eData{idx} = [allECurves(tPoints(idx)-lD,:); allECurves(tPoints(idx) + rD,:)]./N;
    pData{idx} = [allPCurves(tPoints(idx)-lD,:); allPCurves(tPoints(idx) + rD,:)];
end
end
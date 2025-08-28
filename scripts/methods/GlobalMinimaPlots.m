% DESCRIPTION: Generates a figure which illustrates the noise sensitivity of
% global minima for each spiral strip pair. 

% INPUTS: 

% Numspirals - integer specifying the number of spiral strips 

% IntCurves - interaction energies betweeen strips for a given species 

% DVectors - discretised stagger vector

% D - value of D that is used to set the axial periodicity lengthscale,
% that is 5D = g + L. 

function [RobustnessData] = GlobalMinimaPlots(NumSpirals,IntCurves,DVectors,D)
%Load the processed data (these already include the periodic correction and smoothing):
%CollagenData = load('ProcDataRat.mat').OutCell;
%DVectors = CollagenData{1}; IntCurves = CollagenData{2};

%Define pairs of interacting strips: 
SpiralPairs = nchoosek(1:NumSpirals,2);
SpiralPairs = [SpiralPairs;transpose(reshape(repelem(1:NumSpirals,2),[2,NumSpirals]))]; %Include pairs of form n-n

%Obtain the robustness data:
RobustnessData = ComputeMins(DVectors,IntCurves);

%Sort the data in ascending order of energies: 
[RobustnessDSort,SortIdx] = sortrows(RobustnessData,1,'ascend');

%Plot the robustness data: 

%Define the horizontal vector: 
XVec = 1:size(RobustnessData,1);

%Create horizontal axis labels: 
SpiralPairsSorted = SpiralPairs(SortIdx,:);
Names = arrayfun(@(X,Y) strcat(num2str(X),'\textrm{-}',num2str(Y)), SpiralPairsSorted(:,1),SpiralPairsSorted(:,2),'UniformOutput',0);

%Identify minima at 0, D and "other":
MinZero = or(abs(RobustnessDSort(:,2))< 50,abs(RobustnessDSort(:,2)-5*D) < 50);
D_lb = 600; D_ub = 700; %Bounds on what is considered a D-minimum
MinStgrsDiv = arrayfun(@(X) abs(RobustnessDSort(:,2)./X),1:4,'UniformOutput',0);
M = cell2mat(cellfun(@(X) all(and(X>D_lb,X<D_ub),2),MinStgrsDiv,'UniformOutput',0));
MinD = any(M,2); MinOther = ~or(MinD,MinZero);

%Prescribe markers:
Mrkrs = {'ksquare','ko','k^'}; Masks = {MinZero,MinD,MinOther};

%Plot the mean energies with absolute robustness values as error bars:
FontSize = 21; MrkrSize = 10;
fig = figure(4);
T = tiledlayout(2,1);
nexttile(T);
    hold on;
    %Plot different categories of minima: 
    for pltIdx = 1 : numel(Mrkrs)
        plot(XVec(Masks{pltIdx}),RobustnessDSort(Masks{pltIdx},1),...
          Mrkrs{pltIdx},'MarkerFaceColor','k','MarkerSize',MrkrSize)
    end
    %Edit the axes:
    set(gca,'xtick',XVec,'xticklabel',Names,'TickLabelInterpreter','latex',...
        'FontSize',FontSize,'GridAlpha',0.15,'TitleHorizontalAlignment','left');

    %Label the axes:
    title('(A)','Interpreter','latex','FontSize',FontSize)
    ylabel('$\min\limits_{\Delta z}\{\textrm{U}^{p}_{i\textrm{-}j}(\Delta z,g_{\textrm{ref}})\}\,/\, k_{B}T$','Interpreter','latex');
    h(1).TitleHorizontalAlignment = 'left';
    grid("on");

    %Set axis limits:
    ylim([-20,0])
    hold off;

    %Legend: 
    legend('$\Delta z_{\textrm{G}} \approx 0$','$\Delta z_{\textrm{G}} \approx k\textrm{D}$',...
        '$\Delta z_{\textrm{G}} \neq 0,k\textrm{D}$','Interpreter','latex')

%Plot the staggers associated with the mean energy minima:
nexttile(T);
    hold on;
    for pltIdx = 1 : numel(Mrkrs)
        plot(XVec(Masks{pltIdx}),RobustnessDSort(Masks{pltIdx},2),Mrkrs{pltIdx},...
            'MarkerFaceColor','k','MarkerSize',MrkrSize)
    end

    %plot(XVec,RobustnessDSort(:,3),'ko','MarkerFaceColor','k')
    set(gca,'xtick',XVec,'xticklabel',Names,'TickLabelInterpreter','latex',...
        'FontSize',FontSize,'GridAlpha',0.5,'TitleHorizontalAlignment','left');
    ylim([-10,5*D]);
    xlim([1,28]);
    hold on;

    %Label the axes:
    title('(B)','Interpreter','latex','FontSize',FontSize)
    h(2).TitleHorizontalAlignment = 'left';
    ylabel('$\Delta z_{\textrm{G}}$','Interpreter','latex');
    grid on;
    hold off;
    yticks([0,D,2*D,3*D,4*D,5*D]); yticklabels(gca,{'0','\textrm{D}','2\textrm{D}','3\textrm{D}',...
        '4\textrm{D}','5\textrm{D}'})

 %Add global labels: 
 xlabel('Strip pair \,/\, $i\textrm{-}j$','Interpreter','latex','FontSize',FontSize)
end
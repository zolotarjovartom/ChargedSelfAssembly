%DESCRIPTION: Computes the number of charged residues per strip in bovine
%collagen. 
function chargedStats(FPathRes,speciesSplit)
NumStrips = 7;
%Define the paths to the relevant files:
fPath = strcat(FPathRes,{'BulkResData_TypeI|A2|A1|A1_Mammals.mat',...
    'BulkResData_TypeII|A1|A1|A1_Mammals.mat',...
    'BulkResData_TypeIII|A1|A1|A1_Mammals.mat',...
    'BulkResData_TypeV|A1|A2|A1_Mammals.mat',...
    'BulkResData_TypeXI|A3|A2|A1_Mammals.mat',...
    'BulkResData_TypeII|A1|A1|A1_Mammals.mat'});

%Define the plot titles: 
titleStr = {'(A) $\alpha_{2}(\textrm{I})[\alpha_{1}(\textrm{I})]_{2}$',...
    '(B) $[\alpha_{1}(\textrm{II})]_{3}$',...
    '(C) $[\alpha_{1}(\textrm{III})]_{3}$',...
    '(D) $[\alpha_{1}(\textrm{V})]_{2}\alpha_{2}(\textrm{V})$',...
    '(E) $\alpha_{3}(\textrm{XI})\alpha_{2}(\textrm{XI})\alpha_{1}(\textrm{XI})$',...
    strcat('(F)  $[\alpha_{1}(\textrm{II})]_{3}$ \textit{',speciesSplit,'}')};
%Create a plot tiling: 
t = tiledlayout(3,2,'TileIndexing','columnmajor','TileSpacing','Compact');

for tidx = 1 : numel(fPath)
    %Loop over the data:
    nexttile(t);
    bulkData = load(fPath{tidx}).s1;
    
    %Define the set of charged residues: 
    chargedSet = {'2','4','5','7','9','12','19'};
    chargedResInt = [2,4,5,7,9,12,19]';
    
    %Loop over all the data: 
    if tidx < numel(fPath)
        allBarData = zeros(NumStrips,numel(chargedSet),numel(bulkData));
        
        for jidx = 1 : numel(bulkData)
            
            %Extract the residue data per strip: 
            resNameData = bulkData(jidx).StripAAInts;
        
            for idx = 1 : numel(resNameData)
                 %Count the number of each charged residue:
                 allBarData(idx,:,jidx) = histcounts(categorical(cellstr(num2str(resNameData{idx}))),chargedSet);
            end
           
        end
        
        %Calculate the mean residue number: 
        meanResNum = mean(allBarData,3);
        stdResNum = std(allBarData,0,3);
        
        %Plot the bar chart: 
        %X = categorical(strcat({'Strip'},{' '},cellstr(num2str(transpose(1:7)))));
        X = categorical(cellstr(aminolookup(int2aa(chargedResInt))));
        b = bar(X,meanResNum');
        grid on;
        
        %Generate colours for the barchart: 
        cols = distinguishable_colors(NumStrips);
        
        hold on;
        %Amend the bar properties:
        for idx = 1 : NumStrips
            %Bar colour:
            b(idx).FaceColor = cols(idx,:);
        
            %Add errorbars: 
            xEnds = b(idx).XEndPoints;
            errorbar(xEnds,meanResNum(idx,:),stdResNum(idx,:),'LineStyle','none',...
                'Color',[0.3,0.3,0.3],'LineWidth',1)
        end
        hold off;
        
        %Label the axes: 
        ylim([0,35])
        set(gca,'TickLabelInterpreter','latex','FontSize',20,'TitleHorizontalAlignment','left');
    
        %Add titles: 
        title(titleStr{tidx},'Interpreter','latex','Fontsize',20)
        
    %For the last entry, plot just the data for bovine collagen:
    else
        loc = strcmp(bulkData(1).SeqNamesPaired,speciesSplit);
        resNameData = bulkData(loc).StripAAInts;
        specData = zeros(NumStrips,numel(chargedSet));
        for idx = 1 : numel(resNameData)
             %Count the number of each charged residue:
             specData(idx,:) = histcounts(categorical(cellstr(num2str(resNameData{idx}))),chargedSet);
        end
        
        %Plot the data:
         X = categorical(cellstr(aminolookup(int2aa(chargedResInt))));
        b = bar(X,specData');
        grid on;
        
        %Generate colours for the barchart: 
        cols = distinguishable_colors(NumStrips);
        
        hold on;
        %Amend the bar properties:
        for idx = 1 : NumStrips
            %Bar colour:
            b(idx).FaceColor = cols(idx,:);
        end
        hold off;

        %Label the axes: 
        ylim([0,35])
        set(gca,'TickLabelInterpreter','latex','FontSize',20,'TitleHorizontalAlignment','left');
    
        %Add titles: 
        title(titleStr{tidx},'Interpreter','latex','Fontsize',20)
  
    end
end
ylabel(t,'Number of Occurences','Interpreter','latex','FontSize',30);
end
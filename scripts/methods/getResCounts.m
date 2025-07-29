function getResCounts(FPathRes,species,NumStrips)
%Load bulk residue data: 
bulkData = load(FPathRes).s1;

%Define the set of charged residues: 
chargedSet = {'2','4','5','7','9','12','19'};

%Identify the desired species: 
TF_spec = contains(bulkData(1).SeqNamesPaired,species);

%Loop over all the data and count the number of residues: 
allBarData = zeros(NumStrips,numel(chargedSet));
    
for jidx = find(TF_spec)
    
    %Extract the residue data per strip: 
    resNameData = bulkData(jidx).StripAAInts;

    for idx = 1 : numel(resNameData)
         %Count the number of each charged residue:
         allBarData(idx,:) = histcounts(categorical(cellstr(num2str(resNameData{idx}))),chargedSet);
    end
   
end

%Count the net numbers of residues:
resCountsNet = sum(allBarData,1);

%Re-order the residues to be compatible with other scripts: 

%Old order: Arg, Asp, Cys, Glu, His, Lys, Tyr
%New order: Arg, Asp, Glu, Lys, Tyr, Cys, His
%Mapping: 1-->1; 2-->2; 3-->6; 4-->3; 5--> 7; 6-->4; 7-->5
reorderVec = [1,2,4,6,7,3,5];

%Save the residue counts:
resCountsOrdered = resCountsNet(reorderVec);
save('netResCounts|BosTaurus|TypeII|A1|A1|A1.mat','resCountsOrdered')
end
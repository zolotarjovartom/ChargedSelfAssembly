%TITLE: Generate structure 

%DESCRIPTION: given an output from GetBulkResData, produces a structure
%with fieldnames corresponding to the species analysed.

% Generate a structure to store the data: 
function [S,fieldnames] = GenStruct(BulkData)
    %Create a structure which stores the efficient representations:
    RawNames = BulkData.SeqNamesPaired;
    RawNames = regexprep(RawNames,'(\<[a-z])','${upper($1)}');
    fieldnames = cellfun(@(X) erase(X,' '), RawNames,'UniformOutput',0);
    data = zeros(1,size(RawNames,1)); 
    S = cell2struct(num2cell(data,1), fieldnames, 2);
end
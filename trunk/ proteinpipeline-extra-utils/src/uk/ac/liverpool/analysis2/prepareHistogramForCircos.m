%% Creates Histogram files for Circos for mapping peptides. 
%
% Takes as input, a "sorted" file with following fields, where each row 
% lists the locations of an identified peptide -
% Chromosome   start    end     protein-accn  
%   chrIa	17956	18022	gb|TGME49_092920
%   chrIa	18040	18076	gb|TGME49_092920
%   chrIa	18094	18127	gb|TGME49_092920
%   chrIa	18550	19201	gb|TGME49_092920
%   chrIII	1282867	1282903	gb|TGME49_053900
%   chrIII	1282885	1282903	gb|TGME49_053900
%   chrIII	1282921	1282960	gb|TGME49_053900
%   chrIII	1282939	1282960	gb|TGME49_053900
% Example run - prepareHistogramForCircos('chr-Ia.txt',20000)
    
function [records] = prepareHistogramForCircos(chromosomeFile,binSize)

columns = 4;

A = importdata(chromosomeFile);

rows = size(A);

% Find the names of chromosomes in the file
allChrs = {};
for i = 1:rows
    rem = A{i};
    [tok,rem] = strtok(rem);
    allChrs{i} = tok;
end
chrNames = unique(allChrs);

C_records = cell(size(chrNames,2),1);
%for i=1:size(C_records,1)
%    C_records{i} = cell(1,4);
%end

% Read all the data in the file
records = cell(rows,columns);
for i = 1:rows
    rem = A{i};
    for k=1:columns
        [tok,rem] = strtok(rem);
        if(k == 2 || k == 3)
            tok = str2num(tok);
        end
        records{i,k} = tok;
    end
    
    % find index of records{i,1} in C_records and insert the record{i,k}
    % there
    this_chr = records{i,1};
    
    idx = ind2sub(size(chrNames), strmatch(this_chr, chrNames, 'exact'));
    
    content = C_records{idx};
    lenght = size(content,1);
    content{lenght + 1,1} = records{i,1};
    content{lenght + 1,2} = records{i,2};
    content{lenght + 1,3} = records{i,3};
    content{lenght + 1,4} = records{i,4};
    
    C_records{idx} = content;
end

% And then create bin etc for each cell array....

histFile = 'circo-hist-allChr.txt';

for i=1:size(C_records,1)
    chrLabel = chrNames{i};
    bins = createBins(C_records{i},binSize);
    writeToFile(bins,histFile,chrLabel);
end


end


function writeToFile(bins,histFile,chaLabel)

fid = fopen(histFile, 'a');

for i = 1:size(bins,1)
    fprintf(fid,'%s %d %d %d \n',chaLabel ,bins{i,1}, bins{i,2}, bins{i,3});
end
fclose(fid);

end


function [bins] = createBins(records,binSize)

begin = records{1,2};
finish = records{end,3};

% prepare bins
bins = cell(0,3);
for i =1 : binSize: (finish+binSize)
    idx = size(bins,1) + 1;
    bins{idx,1} = i;
    bins{idx,2} = i + binSize - 1;
    bins{idx,3} = 0;
end

% Fill bins
for i=1:size(records,1)
    st = records{i,2};
    ed = records{i,3};
    
    %%%count = records{i,4};
    count  = 1; % we have a pepseq at each line
    
    idx_st = ceil(st/binSize);
    idx_ed = ceil(ed/binSize);
    
    bins{idx_st,3} = bins{idx_st,3} + count;
end

% Do log scale
% for i = 1:size(bins,1)
%     %value = bins{i,3};
%     value = log(bins{i,3});
%     if(isinf(value))
%         value = 0;
%     end
%     bins{i,3} = value;
% end
%     
% Scale counts between 0 and 1
 min_val = min([bins{:,3}]);
 max_val = max([bins{:,3}]);
 range = max_val - min_val;
 
 for i=1:size(bins,1)
     bins{i,3} = (bins{i,3} - min_val) / range;
 end

end


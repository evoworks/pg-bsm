function SEQI = getSequences(data,TAXA,CM,Codon)

lineIDX = strfind(data,newline);
nL = max(setdiff(CM(:,1),CM(:,3)));

SEQI{nL} = [];
discardIDX = [];
for N = 1:nL
    
    idx1 = strfind(data,TAXA{N});
    
    temp = lineIDX - idx1;
    temp(temp < 0) = [];
    idx2 = min(temp) + idx1;
    
    line = data(idx1+length(TAXA{N}):idx2);
    idxT = strfind(line,'T');
    idxC = strfind(line,'C');
    idxA = strfind(line,'A');
    idxG = strfind(line,'G');
    
    line = line(min([idxT,idxC,idxA,idxG]):end);
    line(strfind(line,char(32))) = []; % blanks
    line(strfind(line,char(13))) = []; % carriage return
    line(strfind(line,newline)) = [];  % newlines
    
    n_cod = length(line)/3;
    for site = 1:n_cod
        try
            SEQI{N}(site) = find(strcmp(Codon,line(3*site-2:3*site)) == 1);
        catch
            SEQI{N}(site) = nan;
            discardIDX = [discardIDX;site]; %#ok<AGROW>
        end
    end
    
end

% remove site patterns with missing data

 discardIDX = unique(discardIDX);
 
 for N = 1:nL
     SEQI{N}(discardIDX) = [];
 end

%% END
    
    
    
    
    
    
    
    
    
    
    
    
    

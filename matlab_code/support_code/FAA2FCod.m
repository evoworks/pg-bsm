function Fnew = FAA2FCod(Fold,UAA)

global AminoAcid Codon
num_elements = size(Codon,2);

Fnew = nan*ones(1,num_elements);

for x = 1:20
    
    idx = [];
    
    for y = 1:num_elements
        if strcmp(UAA{x},AminoAcid{y})
            idx = [idx,y]; %#ok<AGROW>
        end
    end
    
   Fnew(idx) = Fold(x);
    
end
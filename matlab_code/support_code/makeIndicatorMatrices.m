function [IS,IN,ID,IID,IIID] = makeIndicatorMatrices(AminoAcid,Codon)

num_elements = size(Codon,2);

IS = zeros(num_elements); % synonymous 
IN = zeros(num_elements); % nonsynonymous
ID = zeros(num_elements); % single nucleotide
IID = zeros(num_elements); % double nucleotide
IIID = zeros(num_elements); % triple nucleotide

% construct indicator matrices

for r = 1:num_elements
    for c = 1:num_elements
        
        if r == c
            continue
        end
        
        dt = sum(diffVector(Codon{r},Codon{c}));
        
        switch dt
            case 1
                ID(r,c) = 1;
            case 2
                IID(r,c) = 1;
            case 3
                IIID(r,c) = 1;
        end
                
        AA1 = findAA(Codon{r},AminoAcid,Codon);
        AA2 = findAA(Codon{c},AminoAcid,Codon);
        
        if strcmp(AA1,AA2) == 1
            IS(r,c) = 1; % synonymous
        else
            IN(r,c) = 1; % non-synonymous
        end
        
    end
end


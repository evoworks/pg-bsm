function pi_nuc = getNucFreq(SEQI,Codon)

num_elements = size(Codon,2);
n_cod = length(SEQI{1});
nL = max(size(SEQI));

% estimate nucleotide frequencies
counts = zeros(3,4);

for n = 1:nL
    for c = 1:n_cod
        
        state = SEQI{n}(c);
         if state > num_elements
            state = state - num_elements;
        end
        
        triplet = Codon{state};
        
        for position = 1:3
            
            counts(position,1) = counts(position,1) + contains(triplet(position),'T');
            counts(position,2) = counts(position,2) + contains(triplet(position),'C');
            counts(position,3) = counts(position,3) + contains(triplet(position),'A');
            counts(position,4) = counts(position,4) + contains(triplet(position),'G');
            
        end
        
    end
end

pi_nuc = counts/(n_cod*nL);





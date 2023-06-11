function TF = makeTargetFrequencies(pi_nuc,Codon)

num_elements = size(Codon,2);

TF = zeros(num_elements);
for r = 1:num_elements
    for c = 1:num_elements
        
        if r == c
            continue
        end
        
        df = diffVector(Codon{r},Codon{c});
        
        tf = 1;
        for position = 1:3
            if df(position) == 1
                idx = strfind('TCAG',Codon{c}(position));
                tf = tf*pi_nuc(position,idx);
            end
        end
        
        TF(r,c) = tf;
        
    end
end


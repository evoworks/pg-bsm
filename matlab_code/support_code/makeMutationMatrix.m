function [M,SNM,DNM,TNM] = makeMutationMatrix(kappa,alpha,beta,TF,Codon)

num_elements = size(TF,1);

M = zeros(num_elements,num_elements);

SNM = zeros(num_elements); % single nucleotide mutation indicator matrix
DNM = zeros(num_elements);
TNM = zeros(num_elements);

for r = 1:num_elements
    for c = 1:num_elements
        
        if r == c
            continue
        end
        
        triplet1 = Codon{r};
        triplet2 = Codon{c};
        
        diff = diffVector(triplet1, triplet2);
        n = sum(diff);
        
        switch n
            case 1
                SNM(r,c) = 1;
            case 2
                DNM(r,c) = 1;
            case 3
                TNM(r,c) = 1;
        end
        
        nt = 0; % number of transitions
        for pos = 1:3
            if strcmpi(triplet1(pos),triplet2(pos))
                continue
            end
            nt = nt + ts_or_tv(triplet1(pos),triplet2(pos));
        end
        
        % construct M(r,c)
        M(r,c) = kappa^nt*TF(r,c);
        
    end
end

M = M.*SNM + alpha*M.*DNM + beta*M.*TNM;

% insert diagonal elements so rows sum to 0
for row = 1:size(M,1)
    M(row,row) = 0;
    M(row,row) = -sum(M(row,:));
end

%% END

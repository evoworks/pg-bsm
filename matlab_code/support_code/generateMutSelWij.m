function [Wij,Sij] = generateMutSelWij(F,Ne)

global Codon
num_elements = size(Codon,2);

% generates a matrix of scaled selection coefficients

Wij = ones(num_elements);
Sij = nan*ones(num_elements);

for i = 1:num_elements
    for j = 1:num_elements
        
        if i == j
            continue
        end
        
        Sij(i,j) = 2*Ne*(F(j) - F(i));
        
        if Sij(i,j) ~= 0
            Wij(i,j) = 2*Sij(i,j)/(1 - exp(-Sij(i,j)));
        end
        
    end
end

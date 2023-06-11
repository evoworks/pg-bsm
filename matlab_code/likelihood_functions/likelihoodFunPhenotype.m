function [LL,L] = likelihoodFunPhenotype(lambda,CM,pi_state,Y)

num_states = length(pi_state);

Q = lambda*ones(num_states)*diag(pi_state);
for row = 1:num_states
    Q(row,row) = 0;
    Q(row,row) = -sum(Q(row,:));
end

Q = Q/dot(pi_state,-diag(Q)); % changed 12 December 2018 

Pt{size(CM,1)} = [];

for branch = 1:size(CM,1)
    Pt{branch} = expm(Q*CM(branch,2));
end

SEQI = cell(length(Y),1);
for n = 1:length(Y)
    SEQI{n} = Y(n);
end

L = PruningAlgorithmW1(pi_state,CM,Pt,SEQI);

LL = -sum(log(L));

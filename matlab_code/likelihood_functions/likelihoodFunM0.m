function [LL,L] = likelihoodFunM0(paramVec,CM,TF,SEQI,AminoAcid,Codon)

% (w, kappa, B)

w = paramVec(1);
kappa = paramVec(2);
B = paramVec(3:end);
CM(:,2) = B(:);

[Q,pi_cod] = makeM0(w,kappa,0,0,TF,AminoAcid,Codon);

% M0 MODEL 

Pt{size(CM,1)} = [];
for branch = 1:size(CM,1)
    Pt{branch} = expm(Q*CM(branch,2));
end

L = PruningAlgorithmW1(pi_cod,CM,Pt,SEQI);

LL = -sum(log(L));


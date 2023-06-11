function [Q,pi_cod] = makeM0(w,kappa,alpha,beta,TF,AminoAcid,Codon)

M = makeMutationMatrix(kappa,alpha,beta,TF,Codon);
P = expm(M*1000); pi_cod = P(1,:);

[IS,IN] = makeIndicatorMatrices(AminoAcid,Codon);
Q = M.*(IS + w*IN);
for row = 1:size(Q,1)
    Q(row,row) = 0;
    Q(row,row) = -sum(Q(row,:));
end

Q = Q/scaleFactor(Q,pi_cod,AminoAcid,Codon);

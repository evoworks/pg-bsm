function [Q1,Q2,pi_cod,r1,r2] = makeM3k2(w1,w2,kappa,alpha,beta,TF,AminoAcid,Codon)

M = makeMutationMatrix(kappa,alpha,beta,TF,Codon);
P = expm(M*1000); pi_cod = P(1,:);

[IS,IN] = makeIndicatorMatrices(AminoAcid,Codon);
Q1 = M.*(IS + w1*IN);
Q2 = M.*(IS + w2*IN);

for row = 1:size(M,1)
    Q1(row,row) = 0; Q1(row,row) = -sum(Q1(row,:));
    Q2(row,row) = 0; Q2(row,row) = -sum(Q2(row,:));
end

r1 = scaleFactor(Q1,pi_cod,AminoAcid,Codon);
r2 = scaleFactor(Q2,pi_cod,AminoAcid,Codon);




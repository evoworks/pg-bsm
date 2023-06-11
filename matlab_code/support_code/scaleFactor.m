function r = scaleFactor(Q,pi_cod,AminoAcid,Codon)

[~,~,ID,IID,IIID] = makeIndicatorMatrices(AminoAcid,Codon);

r1 = sum(sum(diag(pi_cod)*Q.*ID));
r2 = sum(sum(diag(pi_cod)*Q.*IID));
r3 = sum(sum(diag(pi_cod)*Q.*IIID));

r = r1 + 2*r2 + 3*r3;

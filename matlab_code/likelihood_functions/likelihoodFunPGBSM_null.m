function [LL,LPM] = likelihoodFunPGBSM_null(paramVec,CM,TF,SEQI,pi_state,Y,codestr)

eval(['load ' codestr ' AminoAcid Codon'])

% (piM3 w1CL w2CL p1CL delta kappa lambda B) 

piM3 = paramVec(1);
w1M3 = 0.00;
w1CL = paramVec(2);
w2CL = paramVec(3);
p1CL = paramVec(4);
delta = paramVec(5);
kappa = paramVec(6);
lambda = paramVec(7);
B = paramVec(8:end);
CM(:,2) = B(:);

% phenotype component

[~,LPM] = likelihoodFunPhenotype(lambda,CM,pi_state,Y);

% construct rate matrices

[Q1M3,Q1CL,Q2CL,pi_cod,r1M3,r1CL,r2CL] = makeM3k3(w1M3,w1CL,w2CL,kappa,0,0,TF,AminoAcid,Codon);

% scale rate matrices

r = piM3*r1M3 + (1-piM3)*(p1CL*r1CL + (1 - p1CL)*r2CL);

Q1M3 = Q1M3/r;
Q1CL = Q1CL/r;
Q2CL = Q2CL/r;

% M3 likelihood

Pt{size(CM,1)} = [];
for branch = 1:size(CM,1)
    Pt{branch} = expm(Q1M3*CM(branch,2));
end

LM3 = PruningAlgorithmW1(pi_cod,CM,Pt,SEQI);

% CLM3 likelihood

D = diag([p1CL*pi_cod,(1 - p1CL)*pi_cod]);
Z = zeros(size(Codon,2));
I = eye(size(Codon,2));

RCL = [Z (1 - p1CL)*I; p1CL*I Z];
RCL = RCL/sum(sum(D*RCL));

for row = 1:2*size(Codon,2)
    RCL(row,row) = - sum(RCL(row,:));
end

QCL = [Q1CL Z; Z Q2CL] + delta*RCL;

Pt{size(CM,1)} = [];
for branch = 1:size(CM,1)
    Pt{branch} = expm(QCL*CM(branch,2));
end

LCL = PruningAlgorithmW2([p1CL*pi_cod,(1 - p1CL)*pi_cod],CM,Pt,SEQI);

% combined likelihood

L =  [piM3*LM3 + (1-piM3)*LCL,LPM];
LL = -sum(log(L));

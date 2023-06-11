function [LL,L,LM3,LCL] = likelihoodRaMoSS(paramVec,CM,TF,SEQI,codestr,nullModel)

eval(['load ' codestr ' AminoAcid Codon'])

% (propCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta kappa B) 

propCL = paramVec(1);
w1M3 = paramVec(2);
w2M3 = paramVec(3);
p1M3 = paramVec(4);
w1CL = paramVec(5);
w2CL = paramVec(6);
p1CL = paramVec(7);
delta = paramVec(8);
kappa = paramVec(9);
B = paramVec(10:end);
CM(:,2) = B(:);

if nullModel
    delta = 0;
end

% construct rate matrices

% M3 (w1 w2 p1 kappa B)
[Q1M3,Q2M3,pi_cod,r1M3,r2M3] = makeM3k2(w1M3,w2M3,kappa,0,0,TF,AminoAcid,Codon);

% CLM3 (w1 w2 p1 kappa delta B)
[Q1CL,Q2CL,~,r1CL,r2CL] = makeM3k2(w1CL,w2CL,kappa,0,0,TF,AminoAcid,Codon);

% scale rate matrices

r = (1 - propCL)*(p1M3*r1M3 + (1 - p1M3)*r2M3) + propCL*(p1CL*r1CL + (1 - p1CL)*r2CL);

Q1M3 = Q1M3/r;
Q2M3 = Q2M3/r;
Q1CL = Q1CL/r;
Q2CL = Q2CL/r;

% M3 likelihood

Z = zeros(size(Codon,2));
QM3 = [Q1M3 Z; Z Q2M3];

Pt{size(CM,1)} = [];
for branch = 1:size(CM,1)
    Pt{branch} = expm(QM3*CM(branch,2));
end

LM3 = PruningAlgorithmW2([p1M3*pi_cod,(1 - p1M3)*pi_cod],CM,Pt,SEQI);

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

L =  propCL*LCL + (1 - propCL)*LM3;

LL = -sum(log(L));



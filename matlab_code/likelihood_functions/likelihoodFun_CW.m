function [LL,POST,CW_site_likelihoods,LPM]  = likelihoodFun_CW(paramVec,CM,TF,SEQI,pi_state,Y,dMap,dP,codestr)

eval(['load ' codestr ' AminoAcid Codon'])

% (piM3 w1CL w2CL p1CL delta kappa lambda piCW B)

piM3 = paramVec(1);
w1M3 = 0.00;
w1CL = paramVec(2);
w2CL = paramVec(3);
p1CL = paramVec(4);
delta = paramVec(5);
kappa = paramVec(6);
lambda = paramVec(7);
piCW = paramVec(8); % dMap, dP
B = paramVec(9:end);
CM(:,2) = B;

% phenotype component

[~,LPM] = likelihoodFunPhenotype(lambda,CM,pi_state,Y);

% construct rate matrices

[Q1M3,Q1CL,Q2CL,pi_cod,r1M3,r1CL,r2CL] = makeM3k3(w1M3,w1CL,w2CL,kappa,0,0,TF,AminoAcid,Codon);

% scale rate matrices

aaCW = repmat(B(:),1,length(dP));                            % branch lengths
bbCW = repmat(dP(:)',size(dMap,1),1);                        % change map probabilities
ccCW = dMap; ccCW(dMap == 0) = r1CL; ccCW(dMap == 1) = r2CL; % rate ratios

rCW = sum(sum(aaCW.*bbCW.*ccCW))/sum(B);

r = piM3*r1M3 + (1-piM3-piCW)*(p1CL*r1CL + (1 - p1CL)*r2CL) + piCW*rCW;

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

% CW component

sampleSize = length(dP);
n_cod = max(size(SEQI{1}));

Pt{size(CM,1)} = [];
LCW = nan*ones(sampleSize,n_cod);

for sample = 1:sampleSize
    
    for branch = 1:size(CM,1)
        if dMap(branch,sample) == 0
            Pt{branch} = expm(Q1CL*CM(branch,2));
        else
            Pt{branch} = expm(Q2CL*CM(branch,2));
        end
    end
    
    LCW(sample,:) = PruningAlgorithmW1(pi_cod,CM,Pt,SEQI).*(dP(sample)*ones(1,n_cod));
    
end

CW_site_likelihoods = sum(LCW.*repmat(dP(:),1,n_cod));

% combined likelihood

L = piM3*repmat(LM3,sampleSize,1) + (1-piM3-piCW)*repmat(LCL,sampleSize,1) + piCW*LCW;

ell = sum(log(L),2);
max_ell = max(ell);
ell = ell - max_ell;

LL = max_ell + log(dot(dP,exp(ell)));
LL = -sum(LL + log(LPM));


% POST

aa = piM3*LM3;
bb = (1-piM3-piCW)*LCL;
cc = piCW*sum(LCW.*repmat(dP,1,size(LCW,2)),1);
    
POST = [aa(:),bb(:),cc(:)];
POST = POST./repmat(sum(POST,2),1);


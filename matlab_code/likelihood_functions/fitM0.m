function [mle,LL] = fitM0(CM,SEQI,TF,AminoAcid,Codon)

% optimization settings
options = optimoptions('fmincon','algorithm','interior-point','TolX',1e-4,'display','off','UseParallel',true);
options.MaxFunctionEvaluations = 10000;

% vector for balancing root branches
daughter = CM(CM(:,3) == max(CM(:,3)),1);
Aroot = zeros(1,size(CM,1));
Aroot(daughter(1)) = 1;
Aroot(daughter(2)) = -1;

% (w, kappa, B)

disp('Fitting M0...')

lb = [0,1,zeros(1,size(CM,1))];
ub = [50,10,10*ones(1,size(CM,1))];

first_guess = [0.10,4,CM(:,2)'];

Aeq = [0,0,Aroot];
Beq = 0;

tic
f = @(first_guess)likelihoodFunM0(first_guess,CM,TF,SEQI,AminoAcid,Codon);
[mle,LL] = fmincon(f,first_guess,[],[],Aeq,Beq,lb,ub,[],options);
toc



%% END

function [mle,LL] = fitNullPGBSM(CM,SEQI,TF,Y,pi_state,genetic_code)

% optimization settings
options = optimoptions('fmincon','algorithm','interior-point','TolX',1e-4,'display','off','UseParallel',true);
options.MaxFunctionEvaluations = 10000;

% vector for balancing root branches
daughter = CM(CM(:,3) == max(CM(:,3)),1);
Aroot = zeros(1,size(CM,1));
Aroot(daughter(1)) = 1;
Aroot(daughter(2)) = -1;

disp('Fitting null model ...')

% (piM3 w1CL w2CL p1CL delta kappa lambda B)

lb = [0,0,0,0,0,1,0,zeros(1,size(CM,1))];
ub = [1,50,50,1,1,10,1,10*ones(1,size(CM,1))];

first_guess = [0.50,0.10,1.20,0.80,0.20,4,0.2,CM(:,2)'];

A = [0,1,-1,0,0,0,0,zeros(1,size(CM,1))];
B = 0;

Aeq = [0,0,0,0,0,0,0,Aroot];
Beq = 0;

tic
f = @(first_guess)likelihoodFunPGBSM_null(first_guess,CM,TF,SEQI,pi_state,Y,genetic_code);
[mle,LL] = fmincon(f,first_guess,A,B,Aeq,Beq,lb,ub,[],options);
toc

%% END

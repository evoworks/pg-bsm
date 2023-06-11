function [mle0,mle1,LL0,LL1] = fitRaMoSS(MLE_nullPGBM,CM,SEQI,TF,genetic_code)

% optimization settings
options = optimoptions('fmincon','algorithm','interior-point','TolX',1e-4,'display','off','UseParallel',true);
options.MaxFunctionEvaluations = 10000;

% vector for balancing root branches
daughter = CM(CM(:,3) == max(CM(:,3)),1);
Aroot = zeros(1,size(CM,1));
Aroot(daughter(1)) = 1;
Aroot(daughter(2)) = -1;

% fit RaMoSS, mle = (propCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta kappa B)

disp('Fitting RaMoSS with delta = 0 ...')

lb = [0,0,0,0,0,0,0,0,1,zeros(1,size(CM,1))];
ub = [1,50,50,1,50,50,1,1,10,10*ones(1,size(CM,1))];

first_guess = [0.50 0.10 0.50 0.50 0.10 0.50 0.50 0.20 3.00 MLE_nullPGBM(8:end)];

A = [0,1,-1,0,0,0,0,0,0,zeros(1,size(CM,1));
        0,0,0,0,1,-1,0,0,0,zeros(1,size(CM,1))];
B = [0,1];

Aeq = [0,0,0,0,0,0,0,0,0,Aroot];
Beq = 0;

tic
f = @(first_guess)likelihoodRaMoSS(first_guess,CM,TF,SEQI,genetic_code,1);
[mle0,LL0] = fmincon(f,first_guess,A,B,Aeq,Beq,lb,ub,[],options);
toc

mle0(8) = 0; % delta = 0 under the nul model

disp('Fitting RaMoSS with delta estimated...')

lb = [0,0,0,0,0,0,0,0,1,zeros(1,size(CM,1))];
ub = [1,50,50,1,50,50,1,1,10,10*ones(1,size(CM,1))];

first_guess = mle0;
first_guess(8) = 0.20;

tic
f = @(first_guess)likelihoodRaMoSS(first_guess,CM,TF,SEQI,genetic_code,0);
[mle1,LL1] = fmincon(f,first_guess,A,B,Aeq,Beq,lb,ub,[],options);
toc

%% END

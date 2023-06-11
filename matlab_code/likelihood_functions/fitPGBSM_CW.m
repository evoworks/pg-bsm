function [mle,LL,POST] = fitPGBSM_CW(MLE_nullPGBM,CM,SEQI,TF,Y,pi_state,lambda,genetic_code)

% sample the history of the phenotype
[~,~,dMap,dP]= SampleAncestralPhenotypes(lambda,pi_state,CM,Y,1e5);

% optimization settings
options = optimoptions('fmincon','algorithm','interior-point','TolX',1e-4,'display','off','UseParallel',true);
options.MaxFunctionEvaluations = 10000;

% vector for balancing root branches
daughter = CM(CM(:,3) == max(CM(:,3)),1);
Aroot = zeros(1,size(CM,1));
Aroot(daughter(1)) = 1;
Aroot(daughter(2)) = -1;

disp('Fitting the CW model ...')

% (piM3 w1CL w2CL p1CL delta kappa lambda piCW B)
    
lb = [0,0,0,0,0,1,0,0,zeros(1,size(CM,1))];
ub = [1,50,50,1,1,10,1,1,10*ones(1,size(CM,1))];
    
first_guess = [MLE_nullPGBM(1:7),0.10,MLE_nullPGBM(8:end)];
    
A = [0,1,-1,0,0,0,0,0,zeros(1,size(CM,1));
        1,0,0,0,0,0,0,1,zeros(1,size(CM,1))];
B = [0,1];

Aeq = [0,0,0,0,0,0,0,0,Aroot];
Beq = 0;
    
tic
f = @(first_guess)likelihoodFun_CW(first_guess,CM,TF,SEQI,pi_state,Y,dMap,dP,genetic_code);
mle = fmincon(f,first_guess,A,B,Aeq,Beq,lb,ub,[],options);
toc

[LL,POST] = likelihoodFun_CW(mle,CM,TF,SEQI,pi_state,Y,dMap,dP,genetic_code);


%% END

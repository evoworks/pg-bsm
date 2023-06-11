% Generate Alignments using MSmmtDNA
% by C T Jones September 2019

% This code generates alignments using MSmmtDNA (Jones et al. 2018) on a
% user-specified tree with or without changes in landscapes (see setup.m).

% C.T. Jones, N. Youssef, E. Susko, and J.P. Bielawski, Phenomenological
% load on model parameters can lead to false biological conclusions, Mol.
% Biol.Evol. 2018, 35(6):1473-1488

% C.T. Jones, N. Youssef, E. Susko, and J.P. Bielawski, A phenotype-genotype
% codon model for detecting adaptive evolution, Syst. Biol. 2019

clc
close all

global path_to_pgbsm AminoAcid Codon

path_to_pgbsm = pwd;% sets path_to_pgbsm

path_to_data = [path_to_pgbsm '\simulated_alignments\'];

% make paths so Matlab can find things
addpath(path_to_data)
addpath([path_to_pgbsm '\genetic_codes'])
addpath([path_to_pgbsm '\matlab_code\support_code'])

%% load setup data

try
    load setupData % the user must first run simulation_setup.m
catch
    disp('setupData not found: please run simulation_setup.m')
    return
end

%% load mmtDNA data

% load data structures (cell arrays): AminoAcid Codon
load Mammalian_Mitochondrial_GeneticCode

% make a list of unique amino acids
UAA = unique(AminoAcid,'stable');

% determine the number of codon
num_elements = size(Codon,2);

% load mammalian mtDNA data - all of the following were estimated from the
% real alignment of 20 concatenated mammalian mtDNA sequences distributed
% by paml4.9, Yang 2007, PAML4: phylogenetic analysis by maximum likelihood

load MammalianMtData

% SEQI = sequence data, 20 sequences 3331 codons sites each
% TF = F3x4 target frequencies in a 60x60 matrix estimated from SEQI
% kappa = transition/transversion rate ratio estimate from SEQI
% piCOD_Observed = observed codon frequencies in SEQI
% PWF_Observed = observed pairwise codon frequencies in SEQI

%% construct the mutation rate matrix

% constructs the mutation rate matrix - see equation 9 in Jones et al.2018
M = makeMutationMatrix(kappa,alpha,beta,TF);

% computes stationary frequencies under the neutral process
P = expm(M*1000); pi_neutral = P(1,:)/sum(P(1,:));

% the algorithm was tuned to the following values to generate distributions
% of scaled selection coefficients that match those estimated from real
% mammalian mtDNA data (see Jones et al.2018), so the following parameters
% should not be changed

Ne = 1000; u = 0.08; v = 0.02; lb = 0.001; ub = 0.01;

%% draw selection coefficients for the root sequence

[~,~,ID] = makeIndicatorMatrices;

% starting selection coefficients for all sites

sigma_h = nan*ones(1,n_cod);
F = nan*ones(n_cod,num_elements);
for site = 1:n_cod
    
    sigma_h(site) = 0.001 + (0.01-0.001)*betarnd(u,v);
    
    mv = mnrnd(1,piCOD_Observed);
    [~,aa] = ismember(AminoAcid{mv == 1},UAA);
    
    f = zeros(1,20);
    f(aa) = 0.25; idx = setdiff(1:20,aa);
    
    rho = PWF_Observed(aa,idx); rho = rho/max(rho);
    
    f(idx) = rho - abs(mvnrnd(zeros(1,19),eye(19)));
    f = sigma_h(site)*f/std(f);
    
    F(site,:) = FAA2FCod(f,UAA);
    
end

%%  draw selection coefficients

switching_branches = find(changeMap == 1);

if ~isempty(switching_branches)
    
    for b = 1:length(switching_branches)
        
        temp = nan*ones(nBW+nCW+nRW,num_elements);
        
        % BW sites (new f, same sigma)
        for site = 1:nBW
            
            F(site,:) = 0.01*F(site,:)/std(F(site,:));
            
            mv = mnrnd(1,piCOD_Observed);
            [~,aa] = ismember(AminoAcid{mv == 1},UAA);
            
            f = zeros(1,20);
            f(aa) = 0.25; idx = setdiff(1:20,aa);
            
            rho = PWF_Observed(aa,idx); rho = rho/max(rho);
            
            f(idx) = rho - abs(mvnrnd(zeros(1,19),eye(19)));
            
            f = 0.01*f/std(f);
            temp(site,:) = FAA2FCod(f,UAA);
            
        end
        
        % CW sites (same f, smaller sigma)
        for site = nBW+1:nBW+nCW
            
            F(site,:) = 0.01*F(site,:)/std(F(site,:));
            
            f = 0.0005*F(site,:)/std(F(site,:));
            temp(site,:) = FAA2FCod(f,UAA);
            
        end
        
        % RW sites (same f, larger sigma)
        for site = nBW+nCW+1:nBW+nCW+nRW
            
            F(site,:) = 0.0005*F(site,:)/std(F(site,:));
            
            f = 0.01*F(site,:)/std(F(site,:));
            temp(site,:) = FAA2FCod(f,UAA);
            
        end
        
        eval(['F' num2str(switching_branches(b)) '= temp;'])
        
    end
    
end

%% generate rate matrices for sites with fixed landscapes

disp('Generating rate matrices ...')

Sroot = nan*ones(1,n_cod);
A = nan*ones(num_elements,num_elements,n_cod);
r = nan*ones(1,n_cod);

for site = 1:n_cod % fixed sites
    
    A(:,:,site) = M.*generateMutSelWij(F(site,:),Ne);
    
    for row = 1:num_elements
        A(row,row,site) = 0;
        A(row,row,site) = -sum(A(row,:,site));
    end
    
    pi_cod = (1/num_elements)*ones(1,num_elements).*exp(2*Ne*F(site,:));
    pi_cod = pi_cod/sum(pi_cod);
    
    % generate starting codon
    
    idx = find(pi_cod >= 1e-5);
    Sroot(site) = idx(mnrnd(1,pi_cod(idx)/sum(pi_cod(idx))) == 1);
    
    r(site) = scaleFactor(A(:,:,site),pi_cod);
    
end

%% generate rate matrices for sites with peak shifts

if ~isempty(switching_branches)
    
    Aalt = nan*ones(num_elements,num_elements,nBW+nCW+nRW);
    rep = nan*ones(1,nBW+nCW+nRW);
    
    for b = 1:length(switching_branches)
        
        tag = num2str(switching_branches(b));
        
        eval(['Ftag = F' tag ';'])
        
        for site = 1:nBW+nCW+nRW
            
            Aalt(:,:,site) = M.*generateMutSelWij(Ftag(site,:),Ne);
            
            for row = 1:num_elements
                Aalt(row,row,site) = 0;
                Aalt(row,row,site) = -sum(Aalt(row,:,site));
            end
            
            pi_cod = (1/num_elements)*ones(1,num_elements).*exp(2*Ne*Ftag(site,:));
            pi_cod = pi_cod/sum(pi_cod);
            
            rep(site) = scaleFactor(Aalt(:,:,site),pi_cod);
            
        end
        
        eval(['r' tag ' = rep;'])
        eval(['A' tag ' = Aalt;'])
        
    end
    
end

%% compute the scaling factor

temp = nan*r;
for site = 1:n_cod
    
    if and(site >=1, site <= nBW+nCW+nRW)
        
        tot = 0;
        for branch = 1:2*nL-2
            
            if RateMatrixTags(branch,2) == 0
                tot = tot + r(site)*CM(branch,2);
            else
                eval(['tot = tot + r' num2str(RateMatrixTags(branch,2)) '(site)*CM(branch,2);'])
            end
            
        end
        
        temp(site) = tot/sum(CM(:,2));
        
    else
        
        temp(site) = r(site);
        
    end
    
    
end

sf = median(temp);

%% generate alignments

if CM(1,1) == 1
    CM = flipud(CM);
end

for trial = 1:maxTrial
    
    if trial == 1
        
        fid = fopen([path_to_data 'tree_data.txt'],'wt+');
        fprintf(fid,[T ';']);
        fprintf(fid,char(13));
        fclose(fid);
        
        fid = fopen([path_to_data 'phenotype_map.txt'],'wt+');
        fprintf(fid,['phenotypeMap = ' num2str(phenotypeMap') char(13) newline]);
        fprintf(fid,char(13));
        fclose(fid);
        
    end
    
    if trial < 10
        trial_id = ['00' num2str(trial)];
    elseif trial < 100
        trial_id = ['0' num2str(trial)];
    else
        trial_id = num2str(trial);
    end
    
    tic
    disp([num2str(trial) '/' num2str(maxTrial)])
    
    for site = 1:n_cod
        Sroot(site) = SequenceEvolver(Sroot(site),A(:,:,site),50);
    end
    
    % each trial starts at a different root sequence
    eval(['S' num2str(max(CM(:,3))) '= Sroot;'])
    
    for branch = 1:size(CM,1)
        
        eval(['siteIdx_old = S' num2str(CM(branch,3)) ';'])
        siteIdx_new = siteIdx_old;
        
        if ~isempty(switching_branches)
            tag = RateMatrixTags(CM(branch,1),2);
            
            if tag == 0
                
                for n = 1:n_cod %
                    siteIdx_new(n) = SequenceEvolver(siteIdx_old(n),A(:,:,n)/sf,CM(branch,2));
                end
                
            elseif tag ~= 0
                
                eval(['Aepi = A' num2str(tag) ';'])
                
                for n = 1:nBW+nCW+nRW % sites with changed landscapes
                    siteIdx_new(n) = SequenceEvolver(siteIdx_old(n),Aepi(:,:,n)/sf,CM(branch,2));
                end
                
                for n = nBW+nCW+nRW+1:n_cod % sites with static landscapes
                    siteIdx_new(n) = SequenceEvolver(siteIdx_old(n),A(:,:,n)/sf,CM(branch,2));
                end
                
            end
            
        end
        
        eval(['S' num2str(CM(branch,1)) ' = siteIdx_new;'])
        
    end
    
    % construct SEQI
    
    nL = max(setdiff(CM(:,1),CM(:,3))); % number of leaf nodes
    
    SEQI  = [];
    for node = 1:nL
        for n = 1:n_cod
            
            eval(['temp = S' num2str(node) '(' num2str(n) ');'])
            eval(['SEQI{' num2str(node) '}(' num2str(n) ') = temp;'])
            
        end
    end
    
    % generate sequence text file
    
    fid = fopen([path_to_data 'seqfile_' trial_id '.txt'],'wt+');
    fprintf(fid,[num2str(nL) ' ' num2str(3*n_cod) newline char(13)]);
    
    ind = randsample(nL,nL);
    
    for sqs = 1:nL % sequences
        
        fprintf(fid,[Labels{sqs} char(32) char(32)]);
        
        for cds = 1:n_cod % codons
            r = SEQI{sqs}(cds);
            fprintf(fid,[Codon{r} char(32)]);
        end
        
        fprintf(fid,[newline char(13)]);
    end
    
    fclose(fid);
    
    toc
    
end

%% END





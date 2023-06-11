% batch-process alignments simulated using MSmmtDNA
% by C T Jones September 2019

clc
clear
close all

global path_to_pgbsm

path_to_pgbsm = pwd;% sets path_to_pgbsm

path_to_data = [path_to_pgbsm '\simulated_alignments\']; 

% make paths so Matlab can find things
addpath(path_to_data)
addpath([path_to_pgbsm '\genetic_codes'])
addpath([path_to_pgbsm '\matlab_code\support_code'])
addpath([path_to_pgbsm '\matlab_code\likelihood_functions'])

% load Codon and AminoAcid data structures
genetic_code = 'Mammalian_Mitochondrial_GeneticCode';
load(genetic_code);

% get common data
fileList_tree = dir([path_to_data '\tree_data*']);
fileList_phenotypeMap = dir([path_to_data '\phenotype_map*']);
fileList_labels = dir([path_to_data '\taxon_labels*']);

[T,nL,CM] = readTreeFile(path_to_data,fileList_tree.name);
[Y,pi_state] = readPhenotypeMap(path_to_data,fileList_phenotypeMap.name);
Labels = getLabels(path_to_data,fileList_labels.name);

% get a list of files to process
fileList = dir([path_to_data '\seqfile*.txt']); 
Nfiles = size(fileList,1);

%% fit models

Sim_Output(Nfiles).Nul = [];
Sim_Output(Nfiles).BW = [];
Sim_Output(Nfiles).CW = [];
Sim_Output(Nfiles).rCW = [];
Sim_Output(Nfiles).nulRaMoSS = [];
Sim_Output(Nfiles).nulRaMoSS = [];
Sim_Output(Nfiles).altRaMoSS = [];
Sim_Output(Nfiles).altRaMoSS = [];

for file = 1%:Nfiles
    
    % construct alignment data structure
    data = readSequenceFile(path_to_data,fileList(file).name);
    SEQI = getSequences(data,Labels,CM,Codon);
    
    % construct matrix of target frequencies using F3x4
    pi_nuc = getNucFreq(SEQI,Codon);
    TF = makeTargetFrequencies(pi_nuc,Codon);
    
    % fit null model, mle = (piM3 w1CL w2CL p1CL delta kappa lambda B) 
    [null_mle,LL] = fitNullPGBSM(CM,SEQI,TF,Y,pi_state,genetic_code);
    
    Sim_Output(file).Nul.mle = null_mle;
    Sim_Output(file).Nul.LL = -LL;

    lambda = null_mle(5);
    CM(:,2) = null_mle(8:end);
    
    % fit alternative model, branchwise process
    [mle,LL,POST] = fitPGBSM_BW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
    
    Sim_Output(file).BW.mle = mle;
    Sim_Output(file).BW.LL = -LL;
    Sim_Output(file).BW.POST = POST;

    % fit alternative model, cladewise process
    [mle,LL,POST] = fitPGBSM_CW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
    
    Sim_Output(file).CW.mle = mle;
    Sim_Output(file).CW.LL = -LL;
    Sim_Output(file).CW.POST = POST;

    % fit alternative model, reverse cladewise process
    [mle,LL,POST] = fitPGBSM_rCW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
    
    Sim_Output(file).rCW.mle = mle;
    Sim_Output(file).rCW.LL = -LL;
    Sim_Output(file).rCW.POST = POST;
    
    % fit null and alternate RaMoSS
    [mle0,mle1,LL0,LL1] = fitRaMoSS(null_mle,CM,SEQI,TF,genetic_code);
    
    Sim_Output(file).nulRaMoSS.mle = mle0;
    Sim_Output(file).nulRaMoSS.LL = -LL0;
    Sim_Output(file).altRaMoSS.mle = mle1;
    Sim_Output(file).altRaMoSS.LL = -LL1;

    save([path_to_data '\Sim_Output'],'Sim_Output')

end

disp('Done!')

%% END
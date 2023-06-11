% visualize simulation data
% by C T Jones September 2019

clc
clear
close all

global path_to_pgbsm

path_to_pgbsm = pwd;% sets path_to_pgbsm

% make paths so Matlab can find things
addpath([path_to_pgbsm '\genetic_codes'])
addpath([path_to_pgbsm '\matlab_code\support_code'])
addpath([path_to_pgbsm '\matlab_code\likelihood_functions'])

% path to data
path_to_data = [path_to_pgbsm '\simulated_alignments\'];
addpath(path_to_data)

genetic_code = 'Mammalian_Mitochondrial_GeneticCode';
load(genetic_code);

% specify file names
tree = 'tree_data.txt';
labels = 'taxon_labels.txt';
phenotype = 'phenotype_map.txt';

fileList = dir([path_to_data '\seqfile*']);
load([path_to_data '\Sim_Output']) % <- specify the name of the output file

for seq = 1:size(Sim_Output,2)
    
    sequence = fileList(seq).name;
    Output = Sim_Output(seq);
    
    if isempty(Output.Nul)
        continue
    end
    
    idx = strfind(sequence,'_');
    tag = sequence(idx+1:idx+3);
    
    % read in data
    
    [T,nL,CM] = readTreeFile(path_to_data,tree);             % tree structure
    [Y,pi_state] = readPhenotypeMap(path_to_data,phenotype); % phenotype map
    
    TAXA = getLabels(path_to_data,labels);                   % taxa labels
    
    data = readSequenceFile(path_to_data,sequence);          % sequence data
    SEQI = getSequences(data,TAXA,CM,Codon);
    
    pi_nuc = getNucFreq(SEQI,Codon);
    TF = makeTargetFrequencies(pi_nuc,Codon);
    
    % Parameter Set for Each Model (Output.MODEL.mle)
    
    % (1) nul PG-BSM: (pi0 w1 w2 p1 delta kappa lambda B)
    % (2) BW PG-BSM:  (pi0 w1 w2 p1 delta kappa lambda piBW B)
    % (3) CW PG-BSM:  (pi0 w1 w2 p1 delta kappa lambda piCW B)
    % (4) rCW PG-BSM: (pi0 w1 w2 p1 delta kappa lambda pirCW B)
    % (5) nul RaMoSS: (piCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta = 0 kappa B)
    % (6) alt RaMoSS: (piCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta kappa B)
    
    % B = vector of branch length estimates
    
    % log-likelihoods and model contrasts
    
    LL(1) = Output.Nul.LL;
    LL(2) = Output.BW.LL;
    LL(3) = Output.CW.LL;
    LL(4) = Output.rCW.LL;
    LL(5) = Output.nulRaMoSS.LL;
    LL(6) = Output.altRaMoSS.LL;
    
    LLR(1) = 2*(LL(2) - LL(1)); Contrast{1} = '$\pi_{\mbox{BW}} > 0$';
    LLR(2) = 2*(LL(3) - LL(1)); Contrast{2} = '$\pi_{\mbox{CW}} > 0$';
    LLR(3) = 2*(LL(4) - LL(1)); Contrast{3} = '$\pi_{\mbox{rCW}} > 0$';
    LLR(4) = 2*(LL(6) - LL(5)); Contrast{4} = '$\delta > 0$';
    
    LL = LL - max(LL);
    [LL,newOrder] = sort(LL);
    
    X ={'nul PG-BSM','BW','CW','rCW','nul RaMoSS','alt RaMoSS'};
    X = X(newOrder);
    
    figure(1),clf
    bar(LL)
    set(gca,'TickLabelInterpreter','Latex','fontsize',16,'XTickLabels',X)
    title('Log-Likelihood - max(Log-Likeliood)','Interpreter','Latex','fontsize',20)
    ylabel('Log-Likelihood Difference','Interpreter','Latex','fontsize',20)
    set(gca,'Linewidth',1,'XLim',[0.5,6.5])
    grid on,box on,set(gcf,'pos',[320,290,1050,420])
    
    H = figure(1); saveas(H,[path_to_pgbsm '\simulated_alignments\results\fig1_seq' tag '.png'])
    
    figure(2),clf
    bar(LLR)
    set(gca,'TickLabelInterpreter','Latex','fontsize',16,'XTickLabels',Contrast)
    title('Model Contrasts (red line = 1\% decision bound)','Interpreter','Latex','fontsize',20)
    ylabel('Log-Likelihood Ratio','Interpreter','Latex','fontsize',20)
    set(gca,'Linewidth',1,'XLim',[0.5,4.5])
    grid on,box on,set(gcf,'pos',[320,290,1050,420])
    
    hold on
    plot([0.5,4.5],chi2inv(0.99,1)*ones(1,2),'r','linewidth',2)
    hold off
    
    H = figure(2); saveas(H,[path_to_pgbsm '\simulated_alignments\results\fig2_seq' tag '.png'])

    % pos-hoc analysis
    
    figure(3),subplot(311)
    bar(Output.BW.POST(:,3))
    set(gca,'TickLabelInterpreter','Latex','fontsize',16)
    title('Posterior Probability a site is BW','Interpreter','Latex','fontsize',20)
    ylabel('P(BW)','Interpreter','Latex','fontsize',20)
    grid on,box on
    
    subplot(312)
    bar(Output.CW.POST(:,3))
    set(gca,'TickLabelInterpreter','Latex','fontsize',16)
    title('Posterior Probability a site is CW','Interpreter','Latex','fontsize',20)
    ylabel('P(CW)','Interpreter','Latex','fontsize',20)
    grid on,box on
    
    subplot(313)
    bar(Output.rCW.POST(:,3))
    set(gca,'TickLabelInterpreter','Latex','fontsize',16)
    title('Posterior Probability a site is rCW','Interpreter','Latex','fontsize',20)
    ylabel('P(rCW)','Interpreter','Latex','fontsize',20)
    grid on,box on
    
    set(gcf,'pos',[320,290,1050,420])
    
    H = figure(3); saveas(H,[path_to_pgbsm '\simulated_alignments\results\fig3_seq' tag '.png'])

    % branch length comparisons
    
    B = nan*ones(2*nL-2,6);
    legendStr = cell(1,6);
    for n = 1:6
        [B(:,n),legendStr{n}] = getBranchLengths(n,Output);
    end
    
    figure(4)
    plot(B,'o-')
    set(gca,'TickLabelInterpreter','Latex','fontsize',16)
    title('Branch-Length Comparisons','Interpreter','Latex','fontsize',20)
    xlabel('Branch Number','Interpreter','Latex','fontsize',20)
    ylabel('Estimated Length','Interpreter','Latex','fontsize',20)
    set(gca,'XTick',1:(2*nL-2),'Xlim',[1,2*nL-2])
    grid on,box on,set(gcf,'pos',[320,290,1050,420])
    
    LG = legend(legendStr);
    LG.Interpreter = 'Latex';
    LG.Location = 'EastOutside';
    
    H = figure(4); saveas(H,[path_to_pgbsm '\simulated_alignments\results\fig4_seq' tag '.png'])

    % posterior BW sites and branches
    
    [post,I] = sort(Output.BW.POST(:,3),'descend');
    CS = cumsum(1-post);
    
    idx = find(1 - CS < 0);
    idx = idx(1) - 1;
    sites = I(1:idx);
    
    fid = fopen([path_to_data 'results\BW_seq' tag '.txt'],'wt+');
    
    for sqs = 1:nL % sequences
        
        for cds = 1:length(sites)
            r = SEQI{sqs}(sites(cds));
            fprintf(fid,[Codon{r} '(' AminoAcid{r} ') ']);
        end
        
        fprintf(fid,[char(32) TAXA{sqs} '(' num2str(Y(sqs)) ')']);
        fprintf(fid,[char(13),newline]);
        
    end
    
    fclose(fid);
    
    % BW PG-BSM:  (pi0 w1 w2 p1 delta kappa lambda piBW B)
    [B,titleStr] = getBranchLengths(2,Output); CM(:,2) = B;
    
    mle = Output.BW.mle; lambda = 0.10;
    [cMap,cP] = SampleAncestralPhenotypes(lambda,pi_state,CM,Y,1e3);
    
    [cP,I] = sort(cP,'descend');
    cMap = cMap(:,I);
    
    discard_idx = find(cP <= 0.001);
    cP(discard_idx) = [];
    cMap(:,discard_idx) = [];
    
    LLWz = nan*cP;
    for s = 1:length(cP)
        LLWz(s) = -likelihoodFun_BW(mle,CM,TF,SEQI,pi_state,Y,cMap(:,s),1,genetic_code);
    end
    
    Pz = cP(1)/(cP(1) + dot(cP(2:end),exp(LLWz(2:end) - LLWz(1))));
    
    TreeDiagram(CM,cMap(:,1),TAXA,Y,12,1)
    title([titleStr ' Post(z) = ' num2str(round(100*Pz)/100)],'Interpreter','Latex','fontsize',14)
    
    H = figure(5); saveas(H,[path_to_pgbsm '\simulated_alignments\results\fig5_seq' tag '.png'])
    
    % posterior CW sites and branches
    
    [post,I] = sort(Output.CW.POST(:,3),'descend');
    CS = cumsum(1-post);
    
    idx = find(1 - CS < 0);
    idx = idx(1) - 1;
    sites = I(1:idx);
    
    fid = fopen([path_to_data 'results\CW_seq' tag '.txt'],'wt+');
    
    for sqs = 1:nL % sequences
        
        for cds = 1:length(sites)
            r = SEQI{sqs}(sites(cds));
            fprintf(fid,[Codon{r} '(' AminoAcid{r} ') ']);
        end
        
        fprintf(fid,[char(32) TAXA{sqs} '(' num2str(Y(sqs)) ')']);
        fprintf(fid,[char(13),newline]);
        
    end
    
    fclose(fid);
    
    % CW PG-BSM:  (pi0 w1 w2 p1 delta kappa lambda piBW B)
    [B,titleStr] = getBranchLengths(3,Output); CM(:,2) = B;
    
    mle = Output.CW.mle; lambda = 0.10;
    [~,~,cMap,cP] = SampleAncestralPhenotypes(lambda,pi_state,CM,Y,1e3);
    
    [cP,I] = sort(cP,'descend');
    cMap = cMap(:,I);
    
    discard_idx = find(cP <= 0.001);
    cP(discard_idx) = [];
    cMap(:,discard_idx) = [];
    
    LLWz = nan*cP;
    for s = 1:length(cP)
        LLWz(s) = -likelihoodFun_CW(mle,CM,TF,SEQI,pi_state,Y,cMap(:,s),1,genetic_code);
    end
    
    Pz = cP(1)/(cP(1) + dot(cP(2:end),exp(LLWz(2:end) - LLWz(1))));
    
    TreeDiagram(CM,cMap(:,1),TAXA,Y,12,1)
    title([titleStr ' Post(z) = ' num2str(round(100*Pz)/100)],'Interpreter','Latex','fontsize',14)
    
    H = figure(6); saveas(H,[path_to_pgbsm '\simulated_alignments\results\fig6_seq' tag '.png'])
    
    % posterior rCW sites and branches
    
    [post,I] = sort(Output.rCW.POST(:,3),'descend');
    CS = cumsum(1-post);
    
    idx = find(1 - CS < 0);
    idx = idx(1) - 1;
    sites = I(1:idx);
    
    fid = fopen([path_to_data 'results\rCW_seq' tag '.txt'],'wt+');
    
    for sqs = 1:nL % sequences
        
        for cds = 1:length(sites)
            r = SEQI{sqs}(sites(cds));
            fprintf(fid,[Codon{r} '(' AminoAcid{r} ') ']);
        end
        
        fprintf(fid,[char(32) TAXA{sqs} '(' num2str(Y(sqs)) ')']);
        fprintf(fid,[char(13),newline]);
        
    end
    
    fclose(fid);
    
    % rCW PG-BSM:  (pi0 w1 w2 p1 delta kappa lambda piBW B)
    [B,titleStr] = getBranchLengths(4,Output); CM(:,2) = B;
    
    mle = Output.rCW.mle; lambda = 0.10;
    [~,~,cMap,cP] = SampleAncestralPhenotypes(lambda,pi_state,CM,Y,1e3);
    
    [cP,I] = sort(cP,'descend');
    cMap = cMap(:,I);
    
    discard_idx = find(cP <= 0.001);
    cP(discard_idx) = [];
    cMap(:,discard_idx) = [];
    
    LLWz = nan*cP;
    for s = 1:length(cP)
        LLWz(s) = -likelihoodFun_rCW(mle,CM,TF,SEQI,pi_state,Y,cMap(:,s),1,genetic_code);
    end
    
    Pz = cP(1)/(cP(1) + dot(cP(2:end),exp(LLWz(2:end) - LLWz(1))));
    
    TreeDiagram(CM,mod(cMap(:,1)+1,2),TAXA,Y,12,1)
    title([titleStr ' Post(z) = ' num2str(round(100*Pz)/100)],'Interpreter','Latex','fontsize',14)
    
    H = figure(7); saveas(H,[path_to_pgbsm '\simulated_alignments\results\fig7_seq' tag '.png'])
    
    % make flat file
    
    makeFlatFile([path_to_data 'results\'],CM,Output,tag)
    
end

%% END
% Set up a simulation scenario using MSmmtDNA
% by C T Jones September 2019

% The code produces data required to run simulate_alignments.m

% C.T. Jones, N. Youssef, E. Susko, and J.P. Bielawski, Phenomenological
% load on model parameters can lead to false biological conclusions, Mol.
% Biol.Evol. 2018, 35(6):1473-1488

% C.T. Jones, N. Youssef, E. Susko, and J.P. Bielawski, A phenotype-genotype
% codon model for detecting adaptive evolution, Syst. Biol. 2019

clc
clear
close all

global path_to_pgbsm

path_to_pgbsm = pwd;% sets path_to_pgbsm

path_to_data = [path_to_pgbsm '\simulated_alignments\'];
addpath([path_to_pgbsm '\matlab_code\support_code'])

%% BEGIN USER INPUT

% specify your topology and branch lengths

T = '(((1:0.25,2:0.25):0.55,(3:0.30,4:0.25):0.60):1.75,((5:0.25,6:0.25):0.65,(7:0.25,8:0.50):0.50):1.75)';

% indicate the number of taxa

nL = 8;

% make labels for the leaf nodes of the tree

Labels{1} = 'harpsichord';
Labels{2} = 'piano forte';
Labels{3} = 'acoustic guitar';
Labels{4} = 'electric guitar';
Labels{5} = 'trumpet';
Labels{6} = 'trombone';
Labels{7} = 'cornet';
Labels{8} = 'tuba';

% specify branches over which the phenotype changed (0 = no change, 1 = change)

changeMap = zeros(1,2*nL-2);
changeMap(3) = 1;
changeMap(12) = 1;

% specify the number of sites to be generated under each process (see Jones 
% et al. 2019 for a description of the three kinds of generating processes)

nBW = 15; % branchwise sites
nCW = 15; % cladewise sites
nRW = 0;  % reverse cladewise sites

% specify the double-triple mutation rate DT:
% alpha = 0.01795; beta = 0.001355; % 2.8786 D  +  0.1215 T = 3% DT
% alpha = 0.03710; beta = 0.003000; % 5.7689 D  +  0.2610 T = 6% DT

alpha = 0.01795; beta = 0.001355;

% specify the scale of the simulation

maxTrial = 5; % the number of alignment to generate, no more than 999
n_cod = 300;  % the number of codon sites in each alignment

%% END USER INPUT - do not edit anything below this message

% make a text file with taxa names

fid = fopen([path_to_data 'taxon_labels.txt'],'wt+');

for n = 1:nL
    fprintf(fid,Labels{n});
    fprintf(fid,newline);
end

fclose(fid);
        
% this function builds the tree data structure CM
CM = buildTreeData(T,nL);

% this makes rate matrix tags for the alignment generating code 
RateMatrixTags = makeRateMatrixTags(changeMap,CM);

% this makes the phenotype map for the analytic code
phenotypeMap = recodePhenotype(RateMatrixTags(1:nL,2));

% visualize the tree: TreeDiagram(CM,changeMap,Labels,fontsize,offset)
TreeDiagram(CM,changeMap,Labels,phenotypeMap,12,1)

str01 = ['(nBW,nCW,nRW) = (' num2str(nBW) ',' num2str(nCW) ', ' num2str(nRW) ')'];

if and(alpha == 0, beta == 0)
    str02 = 'no DT mutations, ';
elseif and(alpha == 0.01795, beta == 0.001355)
    str02 = '3\% DT mutations, ';
elseif and(alpha == 0.03710, beta == 0.003000)
    str02 = '6\% DT mutations, ';
else
    str02 = 'unknown DT mutation, ';
end
    
str03 = [num2str(maxTrial) ' alignments, '];
str04 = [num2str(n_cod) ' codon sites'];

TL = title({str01,[str02,str03,str04]},'Interpreter','Latex','FontSize',15);
POS = TL.Position;

TL.Position = POS + [0,0.05*POS(2),0];

POS = get(gca,'Position');
set(gca,'Pos',POS - [0,0.50*POS(2),0,0])

% Node labels can sometimes extend beyond the right boundary of the tree 
% diagram. To avoid this you can specify an offset > 0 to extend the right
% boundary. Leaf nodes are numbered 1 to nL from top to bottom. Note that 
% branch number corresponds to the number of the node that the branch leads
% to (i.e., the node farthest to the right).

%% save data

save([path_to_data '\setupData'],'nL','T','CM','Labels','changeMap','RateMatrixTags','phenotypeMap',...
     'maxTrial','n_cod','alpha','beta','nBW','nCW','nRW')
 
disp('Setup Data Saved')


%% END

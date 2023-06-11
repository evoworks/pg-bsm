function Likelihood = PruningAlgorithmW2(pi_state,CM,Pt,SEQI)

%% PRUNING ALGORITHM

nL = max(setdiff(CM(:,1),CM(:,3)));
Likelihood = nan*ones(1,size(SEQI{1},2));

for site = 1:size(SEQI{1},2)
    
    interNode = setdiff(CM(:,3),max(CM(:,3)));
    rootNode = max(CM(:,3));
    
    % module 1: leaf nodes
    
    V = zeros(max(CM(:,3)),size(Pt{1},1));
    
    for node = 1:nL
        
        codon = SEQI{node}(site);
        V(node,codon + [0,size(Pt{1},1)/2]) = 1;
        
    end
    
    % module 2: internal nodes
    
    for n = 1:length(interNode)
        
        daughter = CM(CM(:,3) == interNode(n),1);
        
        P1 = V(daughter(1),:)*Pt{daughter(1)}';
        P2 = V(daughter(2),:)*Pt{daughter(2)}';
        
        V(interNode(n),:) = P1.*P2;
        
    end
    
    % module 3: root node
    
    daughter = CM(CM(:,3) == rootNode,1);
    
    P1 = V(daughter(1),:)*Pt{daughter(1)}';
    P2 = V(daughter(2),:)*Pt{daughter(2)}';
    
    V(rootNode,:) = P1.*P2;
    
    Likelihood(site) = dot(pi_state,V(rootNode,:));
    
end

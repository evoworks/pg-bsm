function [cMap,cP,dMap,dP,rootF] = SampleAncestralPhenotypes(lambda,pi_state,CM,Y,sampleSize)

num_states = length(pi_state);
nL = max(setdiff(CM(:,1),CM(:,3)));

QP = lambda*ones(num_states)*diag(pi_state)/mean(CM(:,2));
for row = 1:num_states
    QP(row,row) = 0;
    QP(row,row) = -sum(QP(row,:));
end

Pt{size(CM,1)} = [];
for branch = 1:size(CM,1)
    Pt{branch} = expm(QP*CM(branch,2));
end

site_pattern = cell(length(Y),1);
for n = 1:length(Y)
    site_pattern{n} = Y(n);
end

V = PruningAlgorithmMatrix(CM,Pt,site_pattern);

rootF = nan*ones(1,sampleSize);
RY = nan*ones(2*nL-1,sampleSize);
for trial = 1:sampleSize
    
    RY(1:nL,trial) = Y;
   
    Pr_root_given_data = (V(2*nL-1,:).*pi_state)/dot(V(2*nL-1,:),pi_state);
    RY(2*nL-1,trial) = find(mnrnd(1,Pr_root_given_data) == 1);
    
    rootF(trial) = RY(2*nL-1,trial);
  
    for branch = 2*nL-2:-1:nL+1
        
        branch_ancestor = RY(CM(branch,3),trial);
        
        Pr_state = (V(branch,:).*Pt{branch}(branch_ancestor,:))/dot(V(branch,:),Pt{branch}(branch_ancestor,:));
        RY(branch,trial) = find(mnrnd(1,Pr_state) == 1);
         
    end

end

cMap = zeros(2*nL-2,sampleSize);

for trial = 1:sampleSize
    for node = 1:2*nL-2
        if RY(node,trial) ~= RY(CM(node,3),trial) % parent
            cMap(node,trial) = 1;
        end
    end
end


C = unique(cMap','rows')';
cvec = zeros(size(C,2),1);

for col = 1:size(C,2)
    for n = 1:sampleSize
        if isequal(C(:,col),cMap(:,n))
            cvec(col) = cvec(col)+1;
        end
    end
end
    
cMap = C; cP = cvec/sum(cvec);

idx = find(cP > 1e-3); % cull low-frequency histories
cMap = cMap(:,idx); 
cP = cP(idx);
cP = cP/sum(cP);

dMap = cMap; dP = cP;
for m = 1:length(cP)
    switching_branches = find(cMap(:,m) == 1);
    Tags = makeTags(CM,switching_branches);
    dMap(Tags(:,2)~=0,m) = 1;
end

%% END


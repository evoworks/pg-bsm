function TreeDiagram(CM,changeMap,Labels,Phenotype,fontsize,offset)

if CM(1,1) ~= 1
    CM = flipud(CM);
end

if isempty(fontsize)
    fontsize = 15;
end

nL = max(setdiff(CM(:,1),CM(:,3)));

if isempty(changeMap)
    changeMap = zeros(1,2*nL-2);
end

X = zeros(2*nL-1,1);
for node = 1:2*nL-2
    X(node,1) = sum(CM(path2root(node,CM,'rooted'),2));
end

Y = zeros(2*nL-1,1);
for node = 1:nL
    Y(node,1) = node;
end

jobComplete = 0;
currentTier = unique(CM(Y ~= 0,3));

while jobComplete == 0
    
    for n = 1:length(currentTier)
        daughters = CM(CM(:,3) == currentTier(n),1);
        
        if prod(Y(daughters,1)) ~= 0
            Y(currentTier(n),1) = mean(Y(daughters,1));
        end
        
    end
    
    if Y(end,1) ~= 0
        jobComplete = 1;
    else
        currentTier = unique(CM(Y ~= 0,3));
    end
    
end


figure(),clf
hold on
for node = 2*nL-1:-1:nL+1
    
    daughters = CM(CM(:,3) == node,1);
    
    plot(X(node),nL-Y(node)+1,'ko','markerfacecolor','k')
    
    plot(X(node)*ones(1,2),[nL-Y(node)+1,nL-Y(daughters(1))+1],'k-','linewidth',1 + changeMap(daughters(1))*2)
    plot(X(node)*ones(1,2),[nL-Y(node)+1,nL-Y(daughters(2))+1],'k-','linewidth',1 + changeMap(daughters(2))*2)
    
    plot([X(node),X(daughters(1))],(nL-Y(daughters(1))+1)*ones(1,2),'k-','linewidth',1 + changeMap(daughters(1))*2)
    plot([X(node),X(daughters(2))],(nL-Y(daughters(2))+1)*ones(1,2),'k-','linewidth',1 + changeMap(daughters(2))*2)
    
end

for node = 1:nL
    plot(X(node),nL-Y(node)+1,'ko','markerfacecolor','k')
end

axis off

for node = 1:2*nL-2
    
    if and(node <= nL, ~isempty(Labels))
        t = text(X(node)+0.05,nL-Y(node)+1,[Labels{node} '(' num2str(Phenotype(node)) ')']);
        t.FontName = 'Times';
        t.FontSize = fontsize;
    else
        t = text(X(node)+0.05,nL-Y(node)+1,num2str(node));
        t.FontName = 'Times';
        t.FontSize = fontsize;
    end
    
end

xlim = get(gca,'XLim');

tree_depth = nan*ones(nL,1);
for n = 1:nL
    tree_depth(n) = sum(CM(path2root(n,CM,'rooted'),2));
end
tree_depth = max(tree_depth);

plot([xlim(1),tree_depth],0.5*ones(1,2),'k-')
plot(xlim(1)*ones(1,2),[0.40,0.60],'k-')
plot(tree_depth*ones(1,2),[0.40,0.60],'k-')

t = text(xlim(1)-0.005,0.01,num2str(0));
t.FontName = 'Times';
t.FontSize = fontsize;

t = text(tree_depth-0.01,0.01,num2str(round(100*tree_depth)/100));
t.FontName = 'Times';
t.FontSize = fontsize;

xlim = get(gca,'XLim');
set(gca,'XLim',[xlim(1),xlim(2) + offset])

hold off

%% END

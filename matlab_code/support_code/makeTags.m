function Tags = makeTags(CM,switching_branches)

Tags = [CM(:,1),zeros(size(CM,1),1)];

if ~isempty(switching_branches)
    
    if CM(1,1) == 1
        CM = flipud(CM);
    end
    
    for branch = 1:size(CM,1)
        
        current_node = CM(branch,1);
        
        % if switch occurs on the branch leading to the current node
        if ~isempty(intersect(current_node,switching_branches))
            
            % identify all nodes that descend from the current node (including itself)
            
            daughter_nodes = [];
            
            for i = 1:size(CM,1)
                
                path2i = path2root(i,CM,'rooted');
                if ~isempty(intersect(current_node,path2i))
                    daughter_nodes = [daughter_nodes,i]; %#ok<AGROW>
                end
                
            end
            
            % update tags for daughter nodes
            
            Tags(daughter_nodes,2) = current_node;
            
        end
        
    end
    
end


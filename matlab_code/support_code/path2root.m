function path2i = path2root(i,CM,type)

n = max(setdiff(CM(:,1),CM(:,3))); % number of leaf nodes


if strcmpi(type,'rooted')
    
    rt_node = 2*n-1;
    
else strcmpi(type,'unrooted')
    
    rt_node = 2*n-2;
    
end

if i == rt_node
    path2i = 0; % at the root
    return
end

idx = find(CM(:,1) == i);

path2i = i;

if CM(idx,3) == rt_node
    at_root = 1;
else
    at_root = 0;
end

while ~at_root
    
    new_i = CM(idx,3);
    path2i = [path2i,new_i]; %#ok<AGROW>
    
    idx = find(CM(:,1) == new_i);
    
    if CM(idx,3) == rt_node
        at_root = 1;
    else
        at_root = 0;
    end
    
end



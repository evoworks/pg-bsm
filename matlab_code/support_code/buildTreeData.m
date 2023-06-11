function CM = buildTreeData(T,nL)

CM = nan*ones(2*nL-2,3); % (node, branch length, ancestor)
doneyet = 0;
currentT = T;

counter = nL;

while ~ doneyet
    
    leftID = strfind(currentT,'(');
    rightID = strfind(currentT,')');
    
    % find left-most bracketed pair (x,y)
    
    idx2 = min(rightID); % first right-hand bracket from left to right
    
    leftID(leftID > idx2) = [];
    idx1 = max(leftID); % nearest left hand bracket to the left
    
    tempstr = currentT(idx1+1:idx2-1);
    commaID = strfind(tempstr,',');
    colonID = strfind(tempstr,':');
    
    daughters(1) = str2double(tempstr(1:colonID(1)-1));
    daughters(2) = str2double(tempstr(commaID(1)+1:colonID(2)-1));
    
    B(1) = str2double(tempstr(colonID(1)+1:commaID(1)-1));
    B(2) = str2double(tempstr(colonID(2)+1:end));
    
    CM(daughters(1),:) = [daughters(1) B(1) counter+1];
    CM(daughters(2),:) = [daughters(2) B(2) counter+1];
      
    segL = currentT(1:idx1-1);
    segR = currentT(idx2+1:end);
    
    newT = [segL,num2str(counter + 1),segR];

    currentT = newT;
    counter = counter + 1;
    
    % disp(currentT)
    
    if ~contains(currentT,'(')
        doneyet = 1;
    end
       
end

%% END






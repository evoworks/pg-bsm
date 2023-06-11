function type = ts_or_tv(n1,n2)

if strcmpi(n1,n2)
    type = nan;
    return
end

type1 = 'purine';
if or(strcmpi(n1,'T'),strcmpi(n1,'C'))
    type1 = 'pyrimidine';
end
   
type2 = 'purine';
if or(strcmpi(n2,'T'),strcmpi(n2,'C'))
    type2 = 'pyrimidine'; 
end
    
type = 0; % transversion
if strcmpi(type1,type2)
    type = 1; % transition
end




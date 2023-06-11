function Snew = SequenceEvolver(Sold,A,b)

Pr = expm(A*b);
m_vec = Pr(Sold,:);

% remove states with small probabilities
idx = find(m_vec >= 1e-10);

% randomly generate the substituting codon
Snew = idx(mnrnd(1,m_vec(idx)/sum(m_vec(idx))) == 1);








function [Y,pi_state,Labels] = recodePhenotype(phenotypeMap)

Y = phenotypeMap;
nL = length(Y);

uP = unique(Y);
num_states = length(uP);

temp = Y;
pi_state = zeros(1,num_states);
for state = 1:num_states
    idx = find(Y == uP(state));
    temp(idx) = state;
    pi_state(state) = length(idx);
end
pi_state = pi_state/length(Y);
Y = temp(:);

Labels{nL} = [];% phenotype labels
for n = 1:nL
    Labels{n} = num2str(Y(n));
end

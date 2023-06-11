function [B,titleStr] = getBranchLengths(Model,Output)

% (1) nul PG-BSM: (pi0 w1 w2 p1 delta kappa lambda B)
% (2) BW PG-BSM: (pi0 w1 w2 p1 delta kappa lambda piBW B)
% (3) CW PG-BSM: (pi0 w1 w2 p1 delta kappa lambda piCW B)
% (4) rCW PG-BSM: (pi0 w1 w2 p1 delta kappa lambda pirCW B)
% (5) nul RaMoSS: (piCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta = 0 kappa B)
% (6) alt RaMoSS: (piCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta kappa B)

switch Model
    
    case 1 % nul PG-BSM
        B = Output.Nul.mle(8:end);
        titleStr = 'null PG-BSM';
    case 2 % BW PG-BSM
        B = Output.BW.mle(9:end);
        titleStr = 'BW PG-BSM';
    case 3 % CW PG-BSM
        B = Output.CW.mle(9:end);
        titleStr = 'CW PG-BSM';
    case 4 % rCW PG-BSM
        B = Output.rCW.mle(9:end);
        titleStr = 'rCW PG-BSM';
    case 5 % nul RaMoSS
        B = Output.nulRaMoSS.mle(10:end);
        titleStr = 'nul RaMoSS';
    case 6 % alt RaMoSS
        B = Output.altRaMoSS.mle(10:end);
        titleStr = 'alt RaMoSS';
        
end

%% END
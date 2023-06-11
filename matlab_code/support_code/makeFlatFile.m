function makeFlatFile(path_to_data,CM,Output,tag)

% (1) nul PG-BSM: (pi0 w1 w2 p1 delta kappa lambda B)
% (2) BW PG-BSM:  (pi0 w1 w2 p1 delta kappa lambda piBW B)
% (3) CW PG-BSM:  (pi0 w1 w2 p1 delta kappa lambda piCW B)
% (4) rCW PG-BSM: (pi0 w1 w2 p1 delta kappa lambda pirCW B)
% (5) nul RaMoSS: (piCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta = 0 kappa B)
% (6) alt RaMoSS: (piCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta kappa B)

%%

nL = max(setdiff(CM(:,1),CM(:,3)));

%%

fid = fopen([path_to_data 'Output_seq' tag '.txt'],'wt+');

fprintf(fid,['PG-BSM Model Fit Output - ' date newline newline]);

fprintf(fid,['Branch Length Estimates' newline newline]);

B = nan*ones(2*nL-2,6);
legendStr = cell(1,6);
for n = 1:6
    [B(:,n),legendStr{n}] = getBranchLengths(n,Output);
end

fprintf(fid,['Daughter    Nul    BW    CW    rCW    nulRaMoSS    altRaMoSS    Parent' newline]);

for row = 1:size(CM,1)
    
    fprintf(fid,[num2str(CM(row,1)) char(9)]);
    
    for m = 1:6
        fprintf(fid,[num2str(round(100*B(row,m))/100) char(9)]);
    end
    
    fprintf(fid,[num2str(CM(row,3)) newline]);
    
end

fprintf(fid,newline);
fprintf(fid,['Log Likelihoods' newline newline]);
fprintf(fid,['Nul: ' num2str(1e4*round(Output.Nul.LL)/1e4) newline]);
fprintf(fid,['BW: ' num2str(1e4*round(Output.BW.LL)/1e4) newline]);
fprintf(fid,['CW: ' num2str(1e4*round(Output.CW.LL)/1e4) newline]);
fprintf(fid,['rCW: ' num2str(1e4*round(Output.rCW.LL)/1e4) newline]);
fprintf(fid,['nulRaMoSS: ' num2str(1e4*round(Output.nulRaMoSS.LL)/1e4) newline]);
fprintf(fid,['altRaMoSS: ' num2str(1e4*round(Output.altRaMoSS.LL)/1e4) newline]);

fprintf(fid,newline);
fprintf(fid,['Nul PG-BSM MLEs' newline]);
fprintf(fid,newline);
fprintf(fid,['pi0     w1     w2     p1     delta     kappa     lambda' newline]);

for m = 1:7
    fprintf(fid,[num2str(round(100*Output.Nul.mle(m))/100) char(9)]);
end

fprintf(fid,newline);
fprintf(fid,newline);
fprintf(fid,['BW PG-BSM MLEs' newline]);
fprintf(fid,newline);
fprintf(fid,['pi0     w1     w2     p1     delta     kappa     lambda     piBW' newline]);

for m = 1:8
    fprintf(fid,[num2str(round(100*Output.BW.mle(m))/100) char(9)]);
end

fprintf(fid,newline);
fprintf(fid,newline);
fprintf(fid,['CW PG-BSM MLEs' newline]);
fprintf(fid,newline);
fprintf(fid,['pi0     w1     w2     p1     delta     kappa     lambda     piCW' newline]);

for m = 1:8
    fprintf(fid,[num2str(round(100*Output.CW.mle(m))/100) char(9)]);
end

fprintf(fid,newline);
fprintf(fid,newline);
fprintf(fid,['rCW PG-BSM MLEs' newline]);
fprintf(fid,newline);
fprintf(fid,['pi0     w1     w2     p1     delta     kappa     lambda     pirCW' newline]);

for m = 1:8
    fprintf(fid,[num2str(round(100*Output.rCW.mle(m))/100) char(9)]);
end

fprintf(fid,newline);
fprintf(fid,newline);
fprintf(fid,['nulRaMoSS MLEs' newline]);
fprintf(fid,newline);
fprintf(fid,['piCL    w1M3    w2M3    p1M3    w1CL    w2CL    p1CL    delta    kappa' newline]);

for m = 1:9
    fprintf(fid,[num2str(round(100*Output.nulRaMoSS.mle(m))/100) char(9)]);
end

fprintf(fid,newline);
fprintf(fid,newline);
fprintf(fid,['altRaMoSS MLEs' newline]);
fprintf(fid,newline);
fprintf(fid,['piCL    w1M3    w2M3    p1M3    w1CL    w2CL    p1CL    delta    kappa' newline]);

for m = 1:9
    fprintf(fid,[num2str(round(100*Output.altRaMoSS.mle(m))/100) char(9)]);
end


fprintf(fid,newline);
fprintf(fid,newline);
fprintf(fid,['BW PG-BSM Posteriors' newline]);
fprintf(fid,newline);
fprintf(fid,['site     P(w=0)        P(w1<->w2)       P(BW)' newline]);

[~,I] =sort(Output.BW.POST(:,3),'descend');

for row = 1:size(Output.BW.POST,1)
    
    fprintf(fid,[num2str(I(row)) char(9)]);
    for n = 1:3
        
       xx = Output.BW.POST(I(row),n);
       
       if xx == 0
        fprintf(fid,['0.000' char(9) char(9)]);
       elseif xx == 1
        fprintf(fid,['1.000' char(9) char(9)]);
       else
        fprintf(fid,[num2str(round(1e3*xx)/1e3) char(9) char(9)]);
       end
       
    end
    
    fprintf(fid,newline);
    
end


fprintf(fid,newline);
fprintf(fid,newline);
fprintf(fid,['CW PG-BSM Posteriors' newline]);
fprintf(fid,newline);
fprintf(fid,['site     P(w=0)        P(w1<->w2)       P(CW)' newline]);

[~,I] =sort(Output.CW.POST(:,3),'descend');

for row = 1:size(Output.CW.POST,1)
    
    fprintf(fid,[num2str(I(row)) char(9)]);
    for n = 1:3
        
       xx = Output.CW.POST(I(row),n);
       
       if xx == 0
        fprintf(fid,['0.000' char(9) char(9)]);
       elseif xx == 1
        fprintf(fid,['1.000' char(9) char(9)]);
       else
        fprintf(fid,[num2str(round(1e3*xx)/1e3) char(9) char(9)]);
       end
       
    end
    
    fprintf(fid,newline);
    
end


fprintf(fid,newline);
fprintf(fid,newline);
fprintf(fid,['rCW PG-BSM Posteriors' newline]);
fprintf(fid,newline);
fprintf(fid,['site     P(w=0)        P(w1<->w2)       P(rCW)' newline]);

[~,I] =sort(Output.rCW.POST(:,3),'descend');

for row = 1:size(Output.rCW.POST,1)
    
    fprintf(fid,[num2str(I(row)) char(9)]);
    for n = 1:3
        
       xx = Output.rCW.POST(I(row),n);
       
       if xx == 0
        fprintf(fid,['0.000' char(9) char(9)]);
       elseif xx == 1
        fprintf(fid,['1.000' char(9) char(9)]);
       else
        fprintf(fid,[num2str(round(1e3*xx)/1e3) char(9) char(9)]);
       end
       
    end
    
    fprintf(fid,newline);
    
end

fclose(fid);

%% END
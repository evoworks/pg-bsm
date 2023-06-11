function [AA,idx] = findAA(triplet,AminoAcid,Codon)

% matches a codon triplet to its amino acid alias

AA = [];

for idx = 1:size(Codon,2)

    if strcmpi(triplet,Codon{idx})
        AA = AminoAcid{idx};
        break
    end
    
end

% stop and start codons are not included in the CODON lookup table

if isempty(AA)
    AA = 'nonsense';
    idx = nan;
end

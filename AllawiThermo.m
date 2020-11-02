function [ deltaH, deltaS, deltaG37 ] = AllawiThermo( oligo )
%Calculates nearest-neighbor model thermodynamics based on Allawi 1997
%Assumes 1M NaCl
%Assumes WC pairing (for now?)

numBases = length(oligo);

deltaH = 0.0;           %kcal/mol
deltaS = -1.4;          %eu + symmetry correction
deltaG37 = 0.4;         %kcal/mol + symmetry correction

%strand initiation
if strcmp(oligo(1), 'G') || strcmp(oligo(1), 'C')
    deltaH = deltaH + 0.1;
    deltaS = deltaS - 2.8;
    deltaG37 = deltaG37 + 0.98;
elseif strcmp(oligo(1), 'A') || strcmp(oligo(1), 'T')
    deltaH = deltaH + 2.3;
    deltaS = deltaS + 4.1;
    deltaG37 = deltaG37 + 1.03;
else
    error('non-canonical terminal base');
end

%strand initiation (terminal)
if strcmp(oligo(numBases), 'G') || strcmp(oligo(numBases), 'C')
    deltaH = deltaH + 0.1;
    deltaS = deltaS - 2.8;
    deltaG37 = deltaG37 + 0.98;
elseif strcmp(oligo(numBases), 'A') || strcmp(oligo(numBases), 'T')
    deltaH = deltaH + 2.3;
    deltaS = deltaS + 4.1;
    deltaG37 = deltaG37 + 1.03;
else
    error('non-canonical terminal base');
end

%nearest-neighbor calculations
for i = 1:numBases-1
   pair = [oligo(i) oligo(i+1)];
   if strcmp(pair, 'AT')
       deltaH = deltaH - 7.2;
       deltaS = deltaS - 20.4;
       deltaG37 = deltaG37 - 0.88;
   elseif strcmp(pair, 'TA')
       deltaH = deltaH - 7.2;
       deltaS = deltaS - 21.3;
       deltaG37 = deltaG37 - 0.58;
   elseif strcmp(pair, 'CG')
       deltaH = deltaH - 10.6;
       deltaS = deltaS - 27.2;
       deltaG37 = deltaG37 - 2.17;
   elseif strcmp(pair, 'GC')
       deltaH = deltaH - 9.8;
       deltaS = deltaS - 24.4;
       deltaG37 = deltaG37 - 2.24;
   elseif strcmp(pair, 'AA') || strcmp(pair, 'TT')
       deltaH = deltaH - 7.9;
       deltaS = deltaS - 22.2;
       deltaG37 = deltaG37 - 1;
   elseif strcmp(pair, 'CA') || strcmp(pair, 'TG')
       deltaH = deltaH - 8.5;
       deltaS = deltaS - 22.7;
       deltaG37 = deltaG37 - 1.45;
   elseif strcmp(pair, 'GT') || strcmp(pair, 'AC')
       deltaH = deltaH - 8.4;
       deltaS = deltaS - 22.4;
       deltaG37 = deltaG37 - 1.44;
   elseif strcmp(pair, 'CT') || strcmp(pair, 'AG')
       deltaH = deltaH - 7.8;
       deltaS = deltaS - 21.0;
       deltaG37 = deltaG37 - 1.28;
   elseif strcmp(pair, 'GA') || strcmp(pair, 'TC')
       deltaH = deltaH - 8.2;
       deltaS = deltaS - 22.2;
       deltaG37 = deltaG37 - 1.3;
   elseif strcmp(pair, 'GG') || strcmp(pair, 'CC')
       deltaH = deltaH - 8.0;
       deltaS = deltaS - 19.9;
       deltaG37 = deltaG37 - 1.84;
   else
       error(['orthogonal base detected at ' i]);
   end
end

end


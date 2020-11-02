
def AllawiThermo(oligo):
    #Calculates nearest-neighbor model thermodynamics based on Allawi 1997
    #Assumes 1M NaCl
    #Assumes WC pairing (for now?)
    numBases = len(oligo)
    deltaH = 0.0           # kcal/mol
    deltaS = -1.4          # eu + symmetry correction
    deltaG37 = 0.4         # kcal/mol + symmetry correction
    #strand initiation
    if oligo[0] == 'G' or oligo[0] == 'C':
        deltaH = deltaH + 0.1
        deltaS = deltaS - 2.8
        deltaG37 = deltaG37 + 0.98
    elif oligo[0] == 'A' or oligo[0] == 'T':
        deltaH = deltaH + 2.3
        deltaS = deltaS + 4.1
        deltaG37 = deltaG37 + 1.03
    else:
        raise ValueError('non-canonical terminal base')
    #strand initiation (terminal)
    if oligo[numBases-1]== 'G' or oligo[numBases-1] == 'C':
        deltaH = deltaH + 0.1;
        deltaS = deltaS - 2.8;
        deltaG37 = deltaG37 + 0.98;
    elif oligo[numBases-1]== 'A' or oligo[numBases-1] == 'T':
        deltaH = deltaH + 2.3;
        deltaS = deltaS + 4.1;
        deltaG37 = deltaG37 + 1.03;
    else:
        raise ValueError('non-canonical terminal base');
    # nearest-neighbor calculations
    for i in range(0,numBases-1):
        pair = oligo[i] + oligo[i+1]
        if pair == 'AT':
            deltaH = deltaH - 7.2
            deltaS = deltaS - 20.4
            deltaG37 = deltaG37 - 0.88
        elif pair == 'TA':
            deltaH = deltaH - 7.2
            deltaS = deltaS - 21.3
            deltaG37 = deltaG37 - 0.58
        elif pair =='CG':
            deltaH = deltaH - 10.6
            deltaS = deltaS - 27.2
            deltaG37 = deltaG37 - 2.17
        elif pair =='GC':
            deltaH = deltaH - 9.8
            deltaS = deltaS - 24.4
            deltaG37 = deltaG37 - 2.24
        elif (pair == 'AA') or (pair == 'TT'):
            deltaH = deltaH - 7.9
            deltaS = deltaS - 22.2
            deltaG37 = deltaG37 - 1
        elif (pair == 'CA') or (pair == 'TG'):
            deltaH = deltaH - 8.5
            deltaS = deltaS - 22.7
            deltaG37 = deltaG37 - 1.45
        elif (pair =='GT') or (pair =='AC'):
            deltaH = deltaH - 8.4
            deltaS = deltaS - 22.4
            deltaG37 = deltaG37 - 1.44
        elif (pair =='CT') or (pair =='AG'):
            deltaH = deltaH - 7.8
            deltaS = deltaS - 21.0
            deltaG37 = deltaG37 - 1.28
        elif (pair =='GA') or (pair =='TC'):
            deltaH = deltaH - 8.2
            deltaS = deltaS - 22.2
            deltaG37 = deltaG37 - 1.3
        elif (pair =='GG') or (pair =='CC'):
            deltaH = deltaH - 8.0
            deltaS = deltaS - 19.9
            deltaG37 = deltaG37 - 1.84
        else:
            raise ValueError('orthogonal base detected at',i)
    return [deltaH, deltaS, deltaG37]

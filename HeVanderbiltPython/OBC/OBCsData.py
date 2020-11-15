import numpy
import pylab

L = numpy.array([10,14,20,24,26,30,32,34,40,42,46,50])
D = numpy.array([0.7495521937409436,0.748133095630968,0.7551880187600327,0.753688557716125,0.690564577980913,0.728280401470924,0.6950535214361949,0.6326486809343935,0.7794025403496738,0.6853505319229404,0.6899417389825953,0.7033161262753732])
F = numpy.array([0.9618361716626087,0.9428258007428983,0.9552063081901584,0.9446269600248529,0.9170289208571015,0.9526852269574598,0.9526852269574598,0.8214302457216206,0.9913968616328988,0.8768807465924278,0.9083384797698022,0.9106027728730962])

L = 1/L[range(6)]
D = D[range(6)]
F = F[range(6)]

#D = D[::-1]
#F = F[::-1]

pylab.figure('OBCs Residue Sums')
pylab.plot(L,F, marker='s', c='green', markeredgecolor='black')
pylab.plot(L,D, marker='o', c='red', markeredgecolor='black')
pylab.plot(L,F-D, marker='o', c='orange', markeredgecolor='black')

pylab.show()

# PBCS

L = [ 6., 10., 14., 18., 22., 26., 30.]
D = [0.60795342, 0.7953829 , 0.7923511 , 0.77841188, 0.76635186, 0.75685108, 0.74937526]
F = [0.62280857, 0.79948367, 0.80042398, 0.88114045, 0.89914115, 0.91187523, 0.96045241]
R = [0.01485515, 0.00410077, 0.00807289, 0.10272856, 0.13278929, 0.15502415, 0.21107715]

pylab.figure('PBCs Drude Weight and Fsum')
pylab.plot(L,F, marker='s', c='green', markeredgecolor='black')
pylab.plot(L,D, marker='o', c='red', markeredgecolor='black')
pylab.plot(L,R, marker='o', c='orange', markeredgecolor='black')

    
pylab.show()
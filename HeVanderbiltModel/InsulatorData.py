import numpy
import pylab

L = numpy.array([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30])

F = numpy.array([0.48883325+0.j, 0.64692946+0.j, 0.81897152+0.j, 0.92080949+0.j,
       0.95949422+0.j, 1.12745578+0.j, 1.17781494+0.j, 1.19502837+0.j,
       1.22225735+0.j, 1.39989434+0.j, 1.41621769+0.j, 1.52926906+0.j,
       1.54998521+0.j, 1.88841546+0.j, 1.92945105+0.j])

L = numpy.array([2,4,6,8,10,12,14,16])

F = numpy.array([0.48883325+0.j, 0.81897152+0.j,
       0.95949422+0.j, 1.17781494+0.j,
       1.22225735+0.j, 1.41621769+0.j,
       1.54998521+0.j, 1.92945105+0.j])

#L = 1/L[range(0,L.shape[0])]

#D = D[::-1]
#F = F[::-1]

pylab.figure('OBCs Residue Sums')
pylab.plot(L,F, marker='s', c='green', markeredgecolor='black')
#pylab.ylim(-0.05, 2.05)
#pylab.tick_params(axis='y', which='both', right=1, left=1, labelleft='on', labelright='on')

#pylab.show()

# PBCS

L = numpy.array([ 3.,  5.,  7.,  9., 11., 13., 15])
R = numpy.array([4.55131625e-04, 1.96451784e-03, 4.95805238e-01, 1.20510579e+00,
       1.42572660e+00, 1.82982516e+00,1.88841546])

#L = 1/L[range(L.shape[0])]

#pylab.figure('PBCs Drude Weight and Fsum')
pylab.plot(L,R, marker='o', c='orange', markeredgecolor='black')
pylab.ylim(-0.05, 2.05)
pylab.tick_params(axis='y', which='both', right=1, left=1, labelleft='on', labelright='on')

    
pylab.show()
import numpy
import pylab

U = numpy.array([0,5,10,15, 20, 25,30,35,40,45,50])
OBCsDw = numpy.array([1.1000000000000019,1.072019286989176,1.0015489704073033, 0.9099205902004113, 0.8134729254089257,0.7215972291863538,0.6384798776882262,0.5652892950786652,0.5017111398541989,0.446819727201596,0.39951711599614426])
NormalizedOBCsDw = OBCsDw / OBCsDw[0]
y = (OBCsDw+NormalizedOBCsDw)/2
error = abs(OBCsDw-NormalizedOBCsDw)
pylab.figure('OBCsVsPBCs@FermiSurface')
pylab.errorbar(U,y,xerr=None, yerr=error,lolims=True, uplims=True,linewidth=1,fmt='-o')
PBCsDw = numpy.array([1.1414982164090361,1.1121608075159524,1.0387518303166179,0.9438219155575537,0.8442291522802757, 0.7494741222636435,0.6637084294537053,0.588041729138189,0.5221119058417523,0.4649605506420024,0.4154685312748977])
NormalizedPBCsDw = PBCsDw / PBCsDw[0]
y = (PBCsDw+NormalizedPBCsDw)/2
error = abs(PBCsDw-NormalizedPBCsDw)
pylab.errorbar(U,y,xerr=None, yerr=error,lolims=True, uplims=True,linewidth=1,fmt='--o',alpha=0.7)
pylab.legend(('OBCs', 'PBCs'))
#pylab.xticks([])  
#pylab.yticks([])  
#pylab.figure('DwVsU0')
#pylab.plot(U,y)
pylab.show()
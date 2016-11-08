
import pylab
import numpy
import math

from numpy import genfromtxt

pylab.close("all")


## data files
#files = [ './output/data/uniform__accumulated_0.txt',
#          './output/data/uniform__accumulated_1.txt',
#          './output/data/0.500000_accumulated_0.txt',
#          './output/data/0.500000_accumulated_1.txt',
#          './output/data/1.000000_accumulated_0.txt',
#          './output/data/1.000000_accumulated_1.txt',
#          './output/data/1.500000_accumulated_0.txt',
#          './output/data/1.500000_accumulated_1.txt',
#          './output/data/2.000000_accumulated_0.txt',
#          './output/data/2.000000_accumulated_1.txt',
#          './output/data/2.500000_accumulated_0.txt',
#          './output/data/2.500000_accumulated_1.txt' ]
#
#fignames = [ './output/figures/uni__0',
#             './output/figures/uni__1',
#             './output/figures/0_5__0',
#             './output/figures/0_5__1',
#             './output/figures/1_0__0',
#             './output/figures/1_0__1',
#             './output/figures/1_5__0',
#             './output/figures/1_5__1',
#             './output/figures/2_0__0',
#             './output/figures/2_0__1',
#             './output/figures/2_5__0',
#             './output/figures/2_5__1' ]

files = [ './output/data/1.800000_accumulated_0.txt',
          './output/data/1.800000_accumulated_1.txt' ]
          
fignames = [ './output/figures/1_8__0',
             './output/figures/1_8__1' ]

for sim in range(len(fignames)):
    pylab.close('all')
    rwdata = genfromtxt(files[sim], delimiter='\t')
    rwdata = rwdata.T
    
    
    ## plot random walks
    pylab.figure(1)
    for rw in rwdata:
        pylab.plot(rw, alpha = 0.8)
    
    rms   = []
    curve = []
    i = 0
    for step in rwdata.T:
        rms.append(math.sqrt(sum(step[:-1]**2) * 1.0 / step.shape[0]))
        curve.append(pow(i, 0.5))
        i = i +1
        
    pylab.plot(rms,   color = 'k', linewidth = 3.0, alpha = 0.6)
    pylab.plot(curve, color = 'k', linewidth = 3.0)
    
    pylab.savefig(fignames[sim] + '.pdf')
    
    
    ## log scale
    pylab.figure(2)
    for rw in rwdata:
        pylab.plot(rw, alpha = 0.8)
    
    pylab.plot(rms,   color = 'k', linewidth = 3.0, alpha = 0.6)
    pylab.plot(curve, color = 'k', linewidth = 3.0) 
    
    pylab.xscale('log')
    pylab.yscale('log')
    
    pylab.savefig(fignames[sim] + '__log.pdf')


## sym log scale
#pylab.figure(3)

#for rw in rwdata:
#    pylab.plot(rw, alpha = 0.8)

#pylab.plot(rms,   color = 'k', linewidth = 3.0, alpha = 0.6)
#pylab.plot(curve, color = 'k', linewidth = 3.0) 

#pylab.xscale('symlog')
#pylab.yscale('symlog')


## plot hist
#linspace = numpy.linspace(-800, 600, 2000)

#pylab.figure(4)
#pylab.hist(rwdata.T[99999,:-1], linspace, log = True, color = 'k', alpha = 0.4)

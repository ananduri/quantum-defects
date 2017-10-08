import numpy as np
from scipy.io import loadmat,savemat
from cat import cattimepar,climb_dCdrpar
from defects import defectstimepar,climb_dDdrpar
import sys
import time
from powerlaw import gpower
from climbdefects import climbR

N=np.int(sys.argv[1])
#N=100
T=np.float(sys.argv[2])
print "T=",T
#s=np.float(sys.argv[3])
s=0

nums=100

numt=10000

t=np.linspace(-T,T,numt)


#d=loadmat('R1_N24_T0.384_s0.0001.mat')
#rspace=d['r'].ravel()


rspace=np.arange(1,5.1,0.1)

start=time.time()
for r in rspace:
	dr=climb_dCdrpar(s,r,N,T,numt)
	print r,dr
end=time.time()
print end-start
#r0=1.0

#start=time.time()
#r=climbR(r0,N,T,numt,s,nums)
#end=time.time()

#print end-start

#g0=gpower(t,T,r[0])
#g1=gpower(t,T,r[-1])

#D0,f0=defectstimepar(g0,N,T,numt)
#D1,f1=defectstimepar(g1,N,T,numt)

#savemat('/scratch/network/ananduri/defects/R1_N{0}_T{1}_s{2}'.format(N,T,s),{'r':r,'N':N,'T':T,'numt':numt,'s':s,'D0':D0,'D1':D1})

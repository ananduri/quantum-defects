import numpy as np
from scipy.io import savemat,loadmat
from climbdefects import climbDefects,climbP
from climbcat import climbCat
from cat import cattimepar
from defects import defectstimepar
import sys
import time
from powerlaw import gpower

N=np.int(sys.argv[1])

#N=np.pi-(np.pi/N)

T=np.float(sys.argv[2])
s=np.float(sys.argv[3])

nums=100

numt=10000

t=np.linspace(-T,T,numt)
#t=np.linspace(0,T,numt)

#initial pulse from PRL
#g=1-(2*t)/T
#g=g*np.arctan(1)
#g0=np.tan(g)

#g0=g0*-500

#t=np.linspace(-T/2,T/2,numt)


r=1.9

g0=gpower(t,T,r)

#g0=np.linspace(2,0,numt)

#dict=loadmat('C29_N100_T25.1_s0.25.mat')
#g0=dict['g']
#g0=g0[-1,:]

start=time.time()
#g=climbDefects(g0,N,T,numt,s,nums)
g=climbCat(g0,N,T,numt,s,nums)
#g=climbFlips(g0,N,T,numt,s,nums)

#g=climbP(g0,N,T/2,numt,s,nums)


end=time.time()

print end-start

C0=cattimepar(g0,N,T,numt)
C1=cattimepar(g[-1,:],N,T,numt)
#D0,f0=defectstimepar(g0,N,T,numt)
#D1,f1=defectstimepar(g[-1,:],N,T,numt)
#F0,f0=flipstimepar(g0,N,T,numt)
#F1,f1=flipstimepar(g[-1,:],N,T,numt)

#D,f,kD=defectstimepar(g0,N,T,numt)


savemat('/scratch/network/ananduri/defects/C1_N{0}_T{1}_s{2}'.format(N,T,s),{'g':g,'N':N,'T':T,'numt':numt,'s':s,'C0':C0,'C1':C1})
#savemat('/scratch/network/ananduri/defects/D__T{1}'.format(N,T,s),{'g':g0,'N':N,'T':T,'numt':numt,'D':D,'kD':kD})

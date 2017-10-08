import numpy as np
from scipy.io import loadmat,savemat
from cat import cattimepar
from defects import defectstimepar,climb_dDdrpar
import sys
import time
from powerlaw import gpower
from climbdefects import climbR
from climbcat import climbRC

N=np.int(sys.argv[1])
T=np.float(sys.argv[2])
s=np.float(sys.argv[3])

nums=100

numt=10000

t=np.linspace(-T,T,numt)

#d=loadmat('R2_N50_T35.4_s0.001.mat')
#r0=d['r'].ravel()
#r0=r0[-1]

r0=5.5

start=time.time()
#r=climbR(r0,N,T,numt,s,nums)
r=climbRC(r0,N,T,numt,s,nums)
end=time.time()

print end-start

g0=gpower(t,T,r[0])
g1=gpower(t,T,r[-1])

#D0,f0,kD0=defectstimepar(g0,N,T,numt)
#D1,f1,kD1=defectstimepar(g1,N,T,numt)
C0=cattimepar(g0,N,T,numt)
C1=cattimepar(g1,N,T,numt)
print C1[-1]
savemat('/scratch/network/ananduri/defects/fidcat/RC1_N{0}_T{1}_s{2}'.format(N,T,s),{'r':r,'N':N,'T':T,'numt':numt,'s':s,'C0':C0,'C1':C1})

import numpy as np

def power(T):
	delta=1./T
	a=np.log(14.7/delta)
	a=delta*a
	a=a/14.7
	a=-np.log(a)
	return a
	
def gpower(t,T,r):
	g=np.abs(t/T)**r
	g=g*np.sign(t)
	g=1-g
	return g

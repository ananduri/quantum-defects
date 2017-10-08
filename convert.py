import numpy as np
import itertools

def K(N):
	kplus=np.linspace((np.pi/N)-np.pi,np.pi-(np.pi/N),num=N)
	kminus=np.linspace(-np.pi,np.pi-(2*np.pi/N),num=N)
	return np.array([kplus,kminus])

def label(m,N): #must have m<=N
	labels={}
	hash=0
	combs=np.array(list(itertools.combinations(np.arange(N),m)))
	size=combs.shape[0]
	for x in np.arange(size):
		labels[hash]=combs[x,:]
		hash+=1
	return labels
	
def mom(ind,m,N):
	ks=K(N)
	if(m%2==0):
		a=0
	if(m%2==1):
		a=1
	labels=label(m,N)
	set=labels[ind]
	k=np.empty([m])
	for x in np.arange(m): 
		k[x]=ks[a,set[x]]
	return k
	
def js(ind,m,N):
	if(m%2==0):
		a=0
	if(m%2==1):
		a=1
	labels=label(m,N)
	set=labels[ind]+1      #+1 or not?
	return set
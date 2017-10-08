import numpy as np
from convert import K
from getuv import getuv,getab,gettheta,getuvlist,getablist,getuvinf,getabinf,getuvinflist,getabinflist
import multiprocessing
from powerlaw import gpower

def cattime(g,N,T,numt):
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,(N/2)],complex)
	v=np.zeros([numt,(N/2)],complex)
	for kind in np.arange(len(kset)):
		uv=getuvinf(kset[kind],g,T,numt)
		u[:,kind]=uv[:,0]
		v[:,kind]=uv[:,1]
	C=np.zeros([numt])
	for q in xrange(numt):
		c=1.
		for kind in xrange(len(kset)):
			c=c*(np.sin(kset[kind]/2)*u[q,kind] + np.cos(kset[kind]/2)*v[q,kind])*(np.sin(kset[kind]/2)*u[q,kind] + np.cos(kset[kind]/2)*v[q,kind]).conj()
		C[q]=np.real(c)
	return C
	
def dCdg(g,N,T,numt):
	# dgds=np.zeros([numt],complex)
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,len(kset)],complex)
	v=np.zeros([numt,len(kset)],complex)
	x=np.zeros([numt,len(kset)],complex)
	y=np.zeros([numt,len(kset)],complex)
	for kind in np.arange(len(kset)):
		uv=getuvinf(kset[kind],g,T,numt)
		u[:,kind]=uv[:,0]
		v[:,kind]=uv[:,1]
		ab=getabinf(kset[kind],g,T,numt)
		x[:,kind]=ab[:,0]
		y[:,kind]=ab[:,1]
	c=1.
	for kind in xrange(len(kset)):
		c=c*(np.sin(kset[kind]/2)*u[-1,kind] + np.cos(kset[kind]/2)*v[-1,kind])*(np.sin(kset[kind]/2)*u[-1,kind] + np.cos(kset[kind]/2)*v[-1,kind]).conj()
	c=np.real(c)
	summand=0.
	for kind in xrange(len(kset)):
		den=(np.sin(kset[kind]/2)*u[-1,kind] + np.cos(kset[kind]/2)*v[-1,kind])*(np.sin(kset[kind]/2)*u[-1,kind] + np.cos(kset[kind]/2)*v[-1,kind]).conj()
		fact1=np.sin(kset[kind]/2)*x[-1,kind]+np.cos(kset[kind]/2)*y[-1,kind]
		fact2=np.sin(kset[kind]/2)*u[-1,kind].conj() + np.cos(kset[kind]/2)*v[-1,kind].conj()
		fact3=y[:,kind].conj()*v[:,kind] - x[:,kind].conj()*u[:,kind]
		num=np.imag(fact1*fact2*fact3)
		summand=summand+(num/den)
	return np.real(-4*c*summand)
	
def gaussmult(dg,T):
	t=np.linspace(-T,T,len(dg))
	f=-(t**2)+(T**2)
	f=f/np.amax(f)
	dg=dg*f
	return dg
	
def cattimepar(g,N,T,numt):
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,(N/2)],complex)
	v=np.zeros([numt,(N/2)],complex)
	
	poolsize=6
	
	masterlist=[]
	for k in kset:
		list=[]
		list.append(k)
		list.append(g)
		list.append(T)
		list.append(numt)
		masterlist.append(list)
	
	pool=multiprocessing.Pool()
	pool_output=pool.map(getuvinflist,masterlist)
	pool.close()
	pool.join()
	
	for z in xrange(len(kset)):
		uv=pool_output[z]
		u[:,z]=uv[:,0]
		v[:,z]=uv[:,1]
	C=np.zeros([numt])
	for q in xrange(numt):
		c=1.
		for kind in xrange(len(kset)):
			c=c*(np.sin(kset[kind]/2)*u[q,kind] + np.cos(kset[kind]/2)*v[q,kind])*(np.sin(kset[kind]/2)*u[q,kind] + np.cos(kset[kind]/2)*v[q,kind]).conj()
		C[q]=np.real(c)
	return C
	
def dCdgpar(g,N,T,numt):
	# dgds=np.zeros([numt],complex)
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,len(kset)],complex)
	v=np.zeros([numt,len(kset)],complex)
	x=np.zeros([numt,len(kset)],complex)
	y=np.zeros([numt,len(kset)],complex)
	
	poolsize=1
	
	masterlist=[]
	for k in kset:
		list=[]
		list.append(k)
		list.append(g)
		list.append(T)
		list.append(numt)
		masterlist.append(list)
	
	pooluv=multiprocessing.Pool()
	pool_outputuv=pooluv.map(getuvinflist,masterlist)
	pooluv.close()
	
	poolab=multiprocessing.Pool()
	pool_outputab=poolab.map(getabinflist,masterlist)
	poolab.close()	
	
	pooluv.join()
	poolab.join()
	
	for z in xrange(len(kset)):
		uv=pool_outputuv[z]
		xy=pool_outputab[z]
		u[:,z]=uv[:,0]
		v[:,z]=uv[:,1]
		x[:,z]=xy[:,0]
		y[:,z]=xy[:,1]
	c=1.
	for kind in xrange(len(kset)):
		c=c*(np.sin(kset[kind]/2)*u[-1,kind] + np.cos(kset[kind]/2)*v[-1,kind])*(np.sin(kset[kind]/2)*u[-1,kind] + np.cos(kset[kind]/2)*v[-1,kind]).conj()
	c=np.real(c)
	summand=0.
	for kind in xrange(len(kset)):
		den=(np.sin(kset[kind]/2)*u[-1,kind] + np.cos(kset[kind]/2)*v[-1,kind])*(np.sin(kset[kind]/2)*u[-1,kind] + np.cos(kset[kind]/2)*v[-1,kind]).conj()
		fact1=np.sin(kset[kind]/2)*x[-1,kind]+np.cos(kset[kind]/2)*y[-1,kind]
		fact2=np.sin(kset[kind]/2)*u[-1,kind].conj() + np.cos(kset[kind]/2)*v[-1,kind].conj()
		fact3=y[:,kind].conj()*v[:,kind] - x[:,kind].conj()*u[:,kind]
		num=np.imag(fact1*fact2*fact3)
		summand=summand+(num/den)
	return np.real(-4*c*summand)
	
def climb_dCdg(s,g,N,T,numt):
	dg=dCdg(g,N,T,numt)
	dg=gaussmult(dg,T)
	return dg
	
def climb_dCdgpar(s,g,N,T,numt):
	dg=dCdgpar(g,N,T,numt)
	dg=gaussmult(dg,T)
	return dg

def climb_dCdrpar(s,r,N,T,numt):
	t=np.linspace(-T,T,numt)
	g=gpower(t,T,r)
	dg=-dCdgpar(g,N,T,numt)
	dr=(np.abs(t/T)**r)*np.sign(t)*np.log(np.abs(t/T)) #took out negative sign, since 1-C is what we want to min
	return np.dot(dg,dr)

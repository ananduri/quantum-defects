import numpy as np
from convert import K
from getuv import getuv,getab,gettheta,getuvlist,getablist,getuvLlist,getabLlist,h1,getuverr,getuverrlist
import multiprocessing
import time
import os
from powerlaw import gpower
import numpy.random

def defects(g,N,T,numt): 
	start=time.time()
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,(N/2)],complex)
	v=np.zeros([numt,(N/2)],complex)
	for kind in np.arange(len(kset)):
		uv=getuv(kset[kind],g,T,numt)
		u[:,kind]=uv[:,0]
		v[:,kind]=uv[:,1]
	term1=2*np.real(u[-1,:].conj()*v[-1,:])
	term1=-term1*np.sin(kset)
	term2=np.absolute(v[-1,:])**2
	term2=-2*np.cos(kset)*term2
	term1=np.sum(term1)
	term2=np.sum(term2)
	numdef=(N/2)+term1+term2
	end=time.time()
	print end-start
	return numdef/N
	




def dPdg(g,k,T,numt):
	dgds=np.zeros([numt],complex)
	#kset=K(N)[0]
	#kset=kset[(N/2):]
	kset=np.array([k])
	u=np.zeros([numt,len(kset)],complex)
	v=np.zeros([numt,len(kset)],complex)
	x=np.zeros([numt,len(kset)],complex)
	y=np.zeros([numt,len(kset)],complex)
	for kind in np.arange(len(kset)):
		uv=getuv(kset[kind],g,T,numt)
		u[:,kind]=uv[:,0]
		v[:,kind]=uv[:,1]
		ab=getab(kset[kind],g,T,numt)
		x[:,kind]=ab[:,0]
		y[:,kind]=ab[:,1]
	term1=u[-1,:].conj()*x[-1,:] + v[-1,:].conj()*y[-1,:]
	term2=u[-1,:].conj()*y[-1,:] + v[-1,:].conj()*x[-1,:]
	term2=-term2*np.sin(kset)
	term3=-2*np.cos(kset)*v[-1,:].conj()*y[-1,:]
	fact1=term1+term2+term3
	for q in np.arange(numt):
		fact2=y[q,:].conj()*v[q,:]-x[q,:].conj()*u[q,:]
		dgds[q]=np.sum(fact1*fact2)
	return -8*np.imag(dgds)





def dDdg(g,N,T,numt):
	dgds=np.zeros([numt],complex)
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,len(kset)],complex)
	v=np.zeros([numt,len(kset)],complex)
	x=np.zeros([numt,len(kset)],complex)
	y=np.zeros([numt,len(kset)],complex)
	for kind in np.arange(len(kset)):
		uv=getuv(kset[kind],g,T,numt)
		u[:,kind]=uv[:,0]
		v[:,kind]=uv[:,1]
		ab=getab(kset[kind],g,T,numt)
		x[:,kind]=ab[:,0]
		y[:,kind]=ab[:,1]
	term1=u[-1,:].conj()*x[-1,:] + v[-1,:].conj()*y[-1,:]
	term2=u[-1,:].conj()*y[-1,:] + v[-1,:].conj()*x[-1,:]
	term2=-term2*np.sin(kset)
	term3=-2*np.cos(kset)*v[-1,:].conj()*y[-1,:]
	fact1=term1+term2+term3
	for q in np.arange(numt):
		fact2=y[q,:].conj()*v[q,:]-x[q,:].conj()*u[q,:]
		dgds[q]=np.sum(fact1*fact2)
	return -8*np.imag(dgds)
	
def defectstime(g,N,T,numt):
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,(N/2)],complex)
	v=np.zeros([numt,(N/2)],complex)
	for kind in np.arange(len(kset)): #this step parallelized
		uv=getuv(kset[kind],g,T,numt)
		u[:,kind]=uv[:,0]
		v[:,kind]=uv[:,1]
	D=np.zeros([numt])
	kD=np.zeros([len(kset),numt],complex)
	for q in np.arange(numt):
		term1=u[q,:].conj()*v[q,:] + v[q,:].conj()*u[q,:]
		term1=-term1*np.sin(kset)
		term2=v[q,:].conj()*v[q,:]
		term2=-2*np.cos(kset)*term2
		kD[:,q]=((1+term1+term2)/N)
		term1=np.sum(term1)
		term2=np.sum(term2)
		numdef=(N/2)+term1+term2
		D[q]=np.real(numdef/N)
	#also calculate fidelity
	f=np.zeros([numt])
	theta=np.zeros([len(kset)])
	for q in np.arange(numt):
		for kind in np.arange(len(kset)):
			theta[kind]=gettheta(kset[kind],g[q])
		fid=u[q,:]*np.cos(theta) + v[q,:]*np.sin(theta)
		f[q]=np.abs(fid.prod())**2
	return D,f,kD

def gaussmult(dg,T):
	# beta=
	# b=T/beta
	t=np.linspace(-T,T,len(dg))
	# dg=dg*np.exp(-(30./(T**2))*((t-(T/2))**2))
	# a=T*((gi/(2.*(beta**2)))+1.-(2./gi))
	# f=((t-T)**2)*np.exp(-((t-a)/b)**2)
	# f=f/np.amax(abs(f))
	# dg=dg*f
	f=-(t**2)+(T**2)
	f=f/np.amax(f)
	dg=dg*f
	return dg


#no gaussmult for now
def climb_dPdg(s,g,k,T,numt):
	dp=dPdg(g,k,T,numt)
	dp=-gaussmult(dp,T)
	return dp


	
def climb_dDdg(s,g,N,T,numt):
	dg=dDdg(g,N,T,numt)
	dg=-gaussmult(dg,T)
	return dg

def climb_dDdgpar(s,g,N,T,numt):
	dg=dDdgpar(g,N,T,numt)
	dg=-gaussmult(dg,T)
	return dg
	
def climb_dDdh(s,h,N,T,numt):
	g=1.-(h**2)
	dh=dDdg(g,N,T,numt)
	dh=-2*h*dh
	dh=-gauss2mult(dh,T)
	return dh

def climb_dDdhpar(s,h,N,T,numt):
	g=1.-(h**2)
	dh=dDdgpar(g,N,T,numt)
	dh=-2*h*dh
	dh=-gauss2mult(dh,T)
	return dh
	
def gauss2mult(dh,T): 
	n=len(dh)
	n=n/2
	# dh1=dh[:n]
	# dh2=dh[n:]
	t1=np.linspace(-T,0,n)
	t2=np.linspace(0,T,n)
	# dh1=dh1*np.exp(-(30./((T/2.)**2))*((t1+(T/2.))**2))
	# dh2=dh2*np.exp(-(30./((T/2.)**2))*((t2-(T/2.))**2))
	f1=((-(t1+(T/2.))**2)+((T/2.)**2))
	f2=((-(t2-(T/2.))**2)+((T/2.)**2))
	f=np.concatenate((f1,f2))
	f=f/np.amax(f)
	dh=dh*f
	return dh
	
def defectstimeh(h,N,T,numt):
	g=1.-(h**2)
	D,f=defectstime(g,N,T,numt)
	return D,f
	
def defectstimehpar(h,N,T,numt):
	g=1.-(h**2)
	D,f=defectstimepar(g,N,T,numt)
	return D,f
	
# def defectspar(g,N,T,numt): 
	# start=time.time()
	# kset=K(N)[0]
	# kset=kset[(N/2):]
	# u=np.zeros([numt,(N/2)],complex)
	# v=np.zeros([numt,(N/2)],complex)
	
	# outq=multiprocessing.Queue()
	
	# def getuvpar(k):
		# return getuv(k,g,T,numt)
	
	# def worker(kind,outq):
		# uv=getuvpar(kset[kind])
		# u=uv[:,0]
		# v=uv[:,1]
		# dict={'kind':kind,'u':u,'v':v}
		# outq.put(dict)
		# if hasattr(os,'getppid'):
			# print 'parent process:',os.getppid()
		# print 'process id:',os.getpid()
		
	# procs=[]
	# for kind in np.arange(len(kset)):
		# p=multiprocessing.Process(target=worker,args=(kind,outq))
		# procs.append(p)
		# p.start()
	
	# for z in np.arange(len(kset)):
		# dict=outq.get()
		# kindex=dict['kind']
		# u[:,kindex]=dict['u']
		# v[:,kindex]=dict['v']
	
	# for p in procs:
		# p.join()
	
	# term1=2*np.real(u[-1,:].conj()*v[-1,:])
	# term1=-term1*np.sin(kset)
	# term2=np.absolute(v[-1,:])**2
	# term2=-2*np.cos(kset)*term2
	# term1=np.sum(term1)
	# term2=np.sum(term2)
	# numdef=(N/2)+term1+term2
	# end=time.time()
	# print end-start
	# return numdef/N
	
def defectstimepar(g,N,T,numt):
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,(N/2)],complex)
	v=np.zeros([numt,(N/2)],complex)
	
	poolsize=8
	
	masterlist=[]
	for k in kset:
		list=[]
		list.append(k)
		list.append(g)
		list.append(T)
		list.append(numt)
		masterlist.append(list)
	
	pool=multiprocessing.Pool(processes=poolsize)
	pool_output=pool.map(getuvlist,masterlist)
	pool.close()
	pool.join()
	
	for z in xrange(len(kset)):
		uv=pool_output[z]
		u[:,z]=uv[:,0]
		v[:,z]=uv[:,1]
	
	D=np.zeros([numt])
	kD=np.zeros([len(kset),numt],complex)
	for q in np.arange(numt):
		term1=u[q,:].conj()*v[q,:] + v[q,:].conj()*u[q,:]
		term1=-term1*np.sin(kset)
		term2=v[q,:].conj()*v[q,:]
		term2=-2*np.cos(kset)*term2
		kD[:,q]=((1+term1+term2)/N)
		term1=np.sum(term1)
		term2=np.sum(term2)
		numdef=(N/2)+term1+term2
		D[q]=np.real(numdef/N)
	#also calculate fidelity
	f=np.zeros([numt])
	theta=np.zeros([len(kset)])
	for q in np.arange(numt):
		for kind in np.arange(len(kset)):
			theta[kind]=gettheta(kset[kind],g[q])
		fid=u[q,:]*np.cos(theta) + v[q,:]*np.sin(theta)
		f[q]=np.abs(fid.prod())**2
	return D,f,kD


def defectstimeparerr(g,N,T,numt,err):
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,(N/2)],complex)
	v=np.zeros([numt,(N/2)],complex)
	
	poolsize=8
	
	masterlist=[]
	for k in kset:
		list=[]
		list.append(k)
		list.append(g)
		list.append(T)
		list.append(numt)
		list.append(err)
		masterlist.append(list)
	
	pool=multiprocessing.Pool(processes=poolsize)
	pool_output=pool.map(getuverrlist,masterlist)
	pool.close()
	pool.join()
	
	for z in xrange(len(kset)):
		uv=pool_output[z]
		u[:,z]=uv[:,0]
		v[:,z]=uv[:,1]
	
	D=np.zeros([numt])
	kD=np.zeros([len(kset),numt],complex)
	for q in np.arange(numt):
		term1=u[q,:].conj()*v[q,:] + v[q,:].conj()*u[q,:]
		term1=-term1*np.sin(kset)
		term2=v[q,:].conj()*v[q,:]
		term2=-2*np.cos(kset)*term2
		kD[:,q]=((1+term1+term2)/N)
		term1=np.sum(term1)
		term2=np.sum(term2)
		numdef=(N/2)+term1+term2
		D[q]=np.real(numdef/N)
	#also calculate fidelity
	f=np.zeros([numt])
	theta=np.zeros([len(kset)])
	for q in np.arange(numt):
		for kind in np.arange(len(kset)):
			theta[kind]=gettheta(kset[kind],g[q])
		fid=u[q,:]*np.cos(theta) + v[q,:]*np.sin(theta)
		f[q]=np.abs(fid.prod())**2
	return D,f,kD


	
def dDdgpar(g,N,T,numt):
	dgds=np.zeros([numt],complex)
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
	
	pooluv=multiprocessing.Pool(processes=poolsize)
	pool_outputuv=pooluv.map(getuvlist,masterlist)
	pooluv.close()
	
	poolab=multiprocessing.Pool()
	pool_outputab=poolab.map(getablist,masterlist)
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
	
	D=np.zeros([numt])
	for q in np.arange(numt):
		term1=u[q,:].conj()*v[q,:] + v[q,:].conj()*u[q,:]
		term1=-term1*np.sin(kset)
		term2=v[q,:].conj()*v[q,:]
		term2=-2*np.cos(kset)*term2
		term1=np.sum(term1)
		term2=np.sum(term2)
		numdef=(N/2)+term1+term2
		D[q]=np.real(numdef/N)
	term1=u[-1,:].conj()*x[-1,:] + v[-1,:].conj()*y[-1,:]
	term2=u[-1,:].conj()*y[-1,:] + v[-1,:].conj()*x[-1,:]
	term2=-term2*np.sin(kset)
	term3=-2*np.cos(kset)*v[-1,:].conj()*y[-1,:]
	fact1=term1+term2+term3
	for q in np.arange(numt):
		fact2=y[q,:].conj()*v[q,:]-x[q,:].conj()*u[q,:]
		dgds[q]=np.sum(fact1*fact2)
	return -8*np.imag(dgds),D

def climb_dDdrpar(s,r,N,T,numt):
	t=np.linspace(-T,T,numt)
	g=gpower(t,T,r)
	dg,D=dDdgpar(g,N,T,numt)
	dg=-dg
	dr=-(np.abs(t/T)**r)*np.sign(t)*np.log(np.abs(t/T))
	return np.dot(dg,dr),D[-1]

def climb_dDdlambpar(s,lamb,g,N,T,numt):
	dg=dDdlambpar(g,lamb,N,T,numt)
	dg=-gaussmult(dg,T)
	return dg

def dDdlambpar(g,lamb,N,T,numt):
	dgds=np.zeros([numt],complex)
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,len(kset)],complex)
	v=np.zeros([numt,len(kset)],complex)
	x=np.zeros([numt,len(kset)],complex)
	y=np.zeros([numt,len(kset)],complex)
	
	poolsize=5
	
	masterlist=[]
	for k in kset:
		list=[]
		list.append(k)
		list.append(g)
		list.append(lamb)
		list.append(N)
		list.append(T)
		list.append(numt)
		masterlist.append(list)
	
	pooluv=multiprocessing.Pool(processes=poolsize)
	pool_outputuv=pooluv.map(getuvLlist,masterlist)
	pooluv.close()
	
	poolab=multiprocessing.Pool()
	pool_outputab=poolab.map(getabLlist,masterlist)
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
	
	term1=u[-1,:].conj()*x[-1,:] + v[-1,:].conj()*y[-1,:]
	term2=u[-1,:].conj()*y[-1,:] + v[-1,:].conj()*x[-1,:]
	term2=-term2*np.sin(kset)
	term3=-2*np.cos(kset)*v[-1,:].conj()*y[-1,:]
	fact1=term1+term2+term3
	for q in np.arange(numt):
		fact2=y[q,:].conj()*u[q,:]-x[q,:].conj()*v[q,:]
		fact2=-4*1j*h1(g[q],N)*np.sin(kset)*fact2
		dgds[q]=np.sum(fact1*fact2)
	return -8*np.imag(dgds)

def defectsLtimepar(g,lamb,N,T,numt):
	kset=K(N)[0]
	kset=kset[(N/2):]
	u=np.zeros([numt,(N/2)],complex)
	v=np.zeros([numt,(N/2)],complex)
	
	poolsize=5
	
	masterlist=[]
	for k in kset:
		list=[]
		list.append(k)
		list.append(g)
		list.append(lamb)
		list.append(T)
		list.append(numt)
		masterlist.append(list)
	
	pool=multiprocessing.Pool(processes=poolsize)
	pool_output=pool.map(getuvLlist,masterlist)
	pool.close()
	pool.join()
	
	for z in xrange(len(kset)):
		uv=pool_output[z]
		u[:,z]=uv[:,0]
		v[:,z]=uv[:,1]
	
	D=np.zeros([numt])
	for q in np.arange(numt):
		term1=u[q,:].conj()*v[q,:] + v[q,:].conj()*u[q,:]
		term1=-term1*np.sin(kset)
		term2=v[q,:].conj()*v[q,:]
		term2=-2*np.cos(kset)*term2
		term1=np.sum(term1)
		term2=np.sum(term2)
		numdef=(N/2)+term1+term2
		D[q]=np.real(numdef/N)
	return D

def noisepulse(g,eta):
	gerr=np.zeros(g.shape)
	for k in xrange(len(g)):
		gerr[k]=g[k] + numpy.random.uniform(-eta/2,eta/2)
	return gerr

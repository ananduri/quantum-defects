import numpy as np
from scipy.integrate import ode
from convert import K

#goes from -T to T

# def getlist(list):
	# k=list[0]
	# g=list[1]
	# T=list[2]
	# numt=list[3]
	# fifth=list[4]
	
def getuverrlist(list):
	k=list[0]
	g=list[1]
	T=list[2]
	numt=list[3]
	err=list[4]
	return getuverr(k,g,T,numt,err)


def getuvlist(list):
	k=list[0]
	g=list[1]
	T=list[2]
	numt=list[3]
	return getuv(k,g,T,numt)
	
def getablist(list):
	k=list[0]
	g=list[1]
	T=list[2]
	numt=list[3]
	return getab(k,g,T,numt)

def getuvLlist(list):
	k=list[0]
	g=list[1]
	lamb=list[2]
	N=list[3]
	T=list[4]
	numt=list[5]
	return getuvL(k,g,lamb,N,T,numt)

def getabLlist(list):
	k=list[0]
	g=list[1]
	lamb=list[2]
	N=list[3]
	T=list[4]
	numt=list[5]
	return getabL(k,g,lamb,N,T,numt)

def getuvinflist(list):
	k=list[0]
	g=list[1]
	T=list[2]
	numt=list[3]
	return getuvinf(k,g,T,numt)
	
def getabinflist(list):
	k=list[0]
	g=list[1]
	T=list[2]
	numt=list[3]
	return getabinf(k,g,T,numt)
	
def getuvinf(k,g,T,numt): 
	dt=(2*T)/float(numt-1)
	init=np.array([0,1]) 
	uvtime=np.zeros([numt,2],complex)
	uvtime[0,:]=init
	l=1
	UV=ode(deriv).set_integrator('zvode',method='adams',rtol=1e-9)
	UV.set_initial_value(init,-T)
	UV.set_f_params(g,k,T,numt)
	while UV.successful() and UV.t < (T):
		UV.integrate(UV.t+dt)
		if l < numt:
			uvtime[l,:]=UV.y
		l+=1
	return uvtime


def getuverr(k,g,T,numt,err): 
	dt=(2*T)/float(numt-1)
	theta=gettheta(k,g[0]+err)
	init=np.array([np.cos(theta),np.sin(theta)]) 
	uvtime=np.zeros([numt,2],complex)
	uvtime[0,:]=init
	l=1
	UV=ode(deriv).set_integrator('zvode',method='adams',rtol=1e-9)
	UV.set_initial_value(init,-T)
	UV.set_f_params(g,k,T,numt)
	while UV.successful() and UV.t < (T):
		UV.integrate(UV.t+dt)
		if l < numt:
			uvtime[l,:]=UV.y
		l+=1
	return uvtime
	


	
def getabinf(k,g,T,numt): 
	dt=(2*T)/float(numt-1)
	init=np.array([-1,0]) 
	uvtime=np.zeros([numt,2],complex)
	uvtime[0,:]=init
	l=1
	UV=ode(deriv).set_integrator('zvode',method='adams',rtol=1e-9)
	UV.set_initial_value(init,-T)
	UV.set_f_params(g,k,T,numt)
	while UV.successful() and UV.t < (T):
		UV.integrate(UV.t+dt)
		if l < numt:
			uvtime[l,:]=UV.y
		l+=1
	return uvtime
	
def getuv(k,g,T,numt): 
	dt=(2*T)/float(numt-1)
	theta=gettheta(k,g[0])
	init=np.array([np.cos(theta),np.sin(theta)]) #change this for getab
	uvtime=np.zeros([numt,2],complex)
	uvtime[0,:]=init
	l=1
	UV=ode(deriv).set_integrator('zvode',method='adams',rtol=1e-15,atol=1e-15)
	UV.set_initial_value(init,-T)
	UV.set_f_params(g,k,T,numt)
	while UV.successful() and UV.t < (T):
		UV.integrate(UV.t+dt)
		if l < numt:
			uvtime[l,:]=UV.y
		l+=1
	return uvtime
		
def getab(k,g,T,numt):
	dt=(2*T)/float(numt-1)
	theta=gettheta(k,g[0])
	init=np.array([-np.sin(theta),np.cos(theta)]) 
	abtime=np.zeros([numt,2],complex)
	abtime[0,:]=init
	l=1
	AB=ode(deriv).set_integrator('zvode',method='adams',rtol=1e-15,atol=1e-15)
	AB.set_initial_value(init,-T)
	AB.set_f_params(g,k,T,numt)
	while AB.successful() and AB.t < (T):
		AB.integrate(AB.t+dt)
		if l < numt:
			abtime[l,:]=AB.y
		l+=1
	return abtime
	
def getuvL(k,g,lamb,N,T,numt): 
	dt=(2*T)/float(numt-1)
	theta=gettheta(k,g[0])
	init=np.array([np.cos(theta),np.sin(theta)]) #change this for getab
	uvtime=np.zeros([numt,2],complex)
	uvtime[0,:]=init
	l=1
	UV=ode(derivL).set_integrator('zvode',method='adams',atol=1e-12)
	UV.set_initial_value(init,-T)
	UV.set_f_params(g,lamb,N,k,T,numt)
	while UV.successful() and UV.t < (T):
		UV.integrate(UV.t+dt)
		if l < numt:
			uvtime[l,:]=UV.y
		l+=1
	return uvtime
		
def getabL(k,g,lamb,N,T,numt):
	dt=(2*T)/float(numt-1)
	theta=gettheta(k,g[0])
	init=np.array([-np.sin(theta),np.cos(theta)]) 
	abtime=np.zeros([numt,2],complex)
	abtime[0,:]=init
	l=1
	AB=ode(derivL).set_integrator('zvode',method='adams',atol=1e-12)
	AB.set_initial_value(init,-T)
	AB.set_f_params(g,lamb,N,k,T,numt)
	while AB.successful() and AB.t < (T):
		AB.integrate(AB.t+dt)
		if l < numt:
			abtime[l,:]=AB.y
		l+=1
	return abtime

def gettheta(k,g):
	arg=np.sin(k)/(g+np.cos(k))
	th1=-.5*np.arctan(arg)
	if np.sign(np.sin(2*th1))==np.sign(np.sin(k)):
		return th1
	else:
		return th1+(np.pi/2)
	
def deriv(t,xy,g,k,T,numt): 
	# map to to index
	frac=np.float(t+T)/(2*T)
	frac=frac*numt
	ind=np.round(frac)
	gmatch=np.array([])
	# print ind
	if (ind<=1):
		g=g[0]
	elif (ind>=numt):
		g=g[-1]
	else:
		g=g[ind-1]
	
	matrix=np.zeros([2,2],complex)
	matrix[0,1]=2j*np.sin(k)
	matrix[1,0]=2j*np.sin(k)
	matrix[1,1]=4j*(np.cos(k)+g)
	# print matrix
	return np.dot(matrix,xy),g
	
def derivL(t,xy,g,lamb,N,k,T,numt): 
	# map to to index
	frac=np.float(t+T)/(2*T)
	frac=frac*numt
	ind=np.round(frac)
	# print ind
	if (ind<=1):
		g=g[0]
		l=lamb[0]
	elif (ind>=numt):
		g=g[-1]
		l=lamb[-1]
	else:
		g=g[ind-1]
		l=lamb[ind-1]
	
	matrix=np.zeros([2,2],complex)
	matrix[0,1]=2j*(np.sin(k)-2*1j*l*h1(g,N)*np.sin(k))
	matrix[1,0]=2j*(np.sin(k)+2*1j*l*h1(g,N)*np.sin(k))
	matrix[1,1]=4j*(np.cos(k)+g)
	# print matrix
	return np.dot(matrix,xy),g

def h1(g,N):
	kset=K(N)[0]
	kset=kset[(N/2):]
	num=np.sin(kset)*np.sin(kset)
	den=((g+np.cos(kset))**2)+(np.sin(kset)*np.sin(kset))
	return 2.*np.sum(num/den)/(4.*N) #2 is due to different time intervals

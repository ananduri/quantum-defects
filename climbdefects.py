import numpy as np
from scipy.integrate import ode
from defects import climb_dDdgpar,climb_dDdrpar,climb_dDdlambpar
from defects import climb_dPdg

def climbDefects(g0,N,T,numt,s,nums):
	ds=s/float(nums-1)
	gtime=np.zeros([nums,numt])
	gtime[0,:]=g0
	l=1
	G=ode(climb_dDdgpar).set_integrator('vode',method='adams')
	G.set_initial_value(g0,0)
	G.set_f_params(N,T,numt)
	while G.successful() and G.t < s:
		G.integrate(G.t+ds)
		if l < nums:
			gtime[l,:]=G.y
		l+=1
	return gtime

def climbR(r0,N,T,numt,s,nums):
	ds=s/float(nums-1)
	rs=np.zeros([nums])
	rs[0]=r0
	l=1
	R=ode(climb_dDdrpar).set_integrator('vode',method='adams')
	R.set_initial_value(r0,0)
	R.set_f_params(N,T,numt)
	while R.successful() and R.t < s:
		R.integrate(R.t+ds)
		if l < nums:
			rs[l]=R.y
			print R.y
		l+=1
	return rs



def climbP(g0,N,T,numt,s,nums):
	ds=s/float(nums-1)
	gtime=np.zeros([nums,numt])
	gtime[0,:]=g0
	l=1
	G=ode(climb_dPdg).set_integrator('vode',method='adams')
	G.set_initial_value(g0,0)
	G.set_f_params(N,T,numt)
	while G.successful() and G.t < s:
		G.integrate(G.t+ds)
		if l < nums:
			gtime[l,:]=G.y
		l+=1
	return gtime








def climbL(lamb0,g,N,T,numt,s,nums):
	ds=s/float(nums-1)
	ltime=np.zeros([nums,numt])
	ltime[0,:]=lamb0
	l=1
	L=ode(climb_dDdlambpar).set_integrator('vode',method='adams')
	L.set_initial_value(lamb0,0)
	L.set_f_params(g,N,T,numt)
	while L.successful() and L.t < s:
		L.integrate(L.t+ds)
		if l < nums:
			ltime[l,:]=L.y
		l+=1
	return ltime 

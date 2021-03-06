import numpy as np
from scipy.integrate import ode
from cat import climb_dCdg,climb_dCdgpar,climb_dCdrpar

def climbCat(g0,N,T,numt,s,nums):
	ds=s/float(nums-1)
	gtime=np.zeros([nums,numt])
	gtime[0,:]=g0
	l=1
	G=ode(climb_dCdgpar).set_integrator('vode',method='adams')
	G.set_initial_value(g0,0)
	G.set_f_params(N,T,numt)
	while G.successful() and G.t < s:
		G.integrate(G.t+ds)
		if l < nums:
			gtime[l,:]=G.y
		l+=1
	return gtime

def climbRC(r0,N,T,numt,s,nums):
	ds=s/float(nums-1)
	rs=np.zeros([nums])
	rs[0]=r0
	l=1
	R=ode(climb_dCdrpar).set_integrator('vode',method='adams')
	R.set_initial_value(r0,0)
	R.set_f_params(N,T,numt)
	while R.successful() and R.t < s:
		R.integrate(R.t+ds)
		if l < nums:
			rs[l]=R.y
			print R.y
		l+=1
	return rs

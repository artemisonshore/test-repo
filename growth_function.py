"""
Nanoom Lee (nanoom.lee@stonybrook.edu)
This is a part of codes for my research (Cosmology).
Given cosmological parameters, this code calculates the matter growth function solving 2nd order differential equation.
nanoom.py is for numerical derivative and integration.
If have any question, feel free to contact me.
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import scipy.special
import file_open as fo
import nanoom as nn

def run_D (m, w):

	# Parameters
	m_s = (299792458)**-1	# m/s = 1/(3x10^8)
	l = 100
	c = 1
	h = 0.7
	omega_l = 0.73
	eV_to_K = 11604.5250061657
	m_v = m * eV_to_K
	H_0 = 10**5 * m_s * h
	rho_cr = 8.089*10**-11 * h**2 * eV_to_K**4	# K^4
	a_0 = 10**-1
	n = 1000

	scale_factor = np.logspace (np.log10(a_0), 0, n)
	scale_factor = np.array (scale_factor)
	scale_factor_reverse = scale_factor[::-1]
	redshift = 1/scale_factor_reverse - 1

	# Calculate rho_cr using rho_0
	T = 1.95
	E_tot = integrate.quad (lambda p: 2*p**2*np.sqrt(p**2 + m_v**2)/(np.e**(p/T)+1) / (2*np.pi**2) , 0, T*708)
	rho_0 = E_tot[0]
	omega_v = rho_0 / rho_cr
	omega_m = 0.27

	D = []
	D.append (a_0)
	u = []
	u.append (1)
	
	rho_0_list = []
	x = np.arange (10**-1, 1, 10**-2)
	for i in x:
		T = 1.95 / i
		E = integrate.quad(lambda p: 2*p**2*np.sqrt(p**2 + m_v**2)/(np.e**(p/T)+1) / (2*np.pi**2), 0, T*708)
		rho_0_list.append (E[0])

	f = np.interp (np.log(scale_factor), np.log(x), np.log(rho_0_list))
	rho_interpolated = np.e**f 
	rho_interpolated_reverse = rho_interpolated [::-1]

	drho_reverse = nn.derivative (scale_factor_reverse, rho_interpolated_reverse)
	drho = nn.derivative (scale_factor, rho_interpolated)

	if m == 0:
		omega_v = 0
		omega_m = 0.27
		rho_interpolated *= 0
		rho_interpolated_reverse *= 0
		drho_reverse *= 0
		drho *= 0
	
	E = omega_m / (scale_factor)**3 + omega_l / (scale_factor)**(3*(1+w)) + rho_interpolated / rho_cr
	diff_logE = (-3*omega_m / (scale_factor)**4 - 3*(1+w)* omega_l * (scale_factor)**(-4-3*w) + drho/rho_cr) / E

	# Numerical Calculation for Growth Function
	for i in np.arange (0, len(scale_factor)-1):
		a = scale_factor[i]
		da = scale_factor[i+1] - scale_factor[i]
		u_prime = - (3/a + 1/2*diff_logE[i]) * u[i] + 3*(omega_m/a**3 + rho_interpolated[i]/rho_cr) / ( 2*a**2*E[i] ) * D[i]
		D.append ( D[i] + u[i]*da )
		u.append ( u[i] + u_prime*da )
	
	D = np.array (D)
	
	# D(z), d[D(z)(1+z)]/dz, H(z), chi(z)
	D = D / D[-1]
	D_reverse = D[::-1]
	D_Z = D_reverse * (redshift + 1)
	diff_D_Z = nn.derivative (redshift, D_Z)
	hubble = H_0 * np.sqrt (omega_m * (1+redshift)**3 + omega_l * (redshift+1)**(3*(1+w)) + rho_interpolated_reverse / rho_cr)
	chi = []
	for i in range (len(hubble)):
		chi.append (nn.integrate (redshift[:i+1], 1/hubble))
	chi = np.array (chi)

	return h, H_0, omega_m, w, c, m, redshift, D, D_reverse, diff_D_Z, chi, hubble


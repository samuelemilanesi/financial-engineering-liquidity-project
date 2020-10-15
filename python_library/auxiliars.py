"""------------------------------------------------------------
# LIBRARY OF AUXILIAR FUNCTIONS:
#	- normal_vols_swaption_MHJM
#	- functions for Jamshidian trick
#	- buildStruct
#	- create_vector_payment_dates
#	- compute_pi_ul
#	- bootstrap_Z_curve
# 
------------------------------------------------------------"""


#------------------------------------------------------------
#	Used libraries
#------------------------------------------------------------
import pandas as pd
import numpy as np
from datetime import datetime

import utility as ut
import dates_manipulation as dm

#------------------------------------------------------------
#	normal_vols_swaption_MHJM 
# returns reciver cash-settle swaption market prices 
#------------------------------------------------------------

def normal_vols_swaption_MHJM(dates_semiannual,deltas_fixed_leg, discount_factors_semiannual,
                              forward_beta_semiannual, settlement_swaption, nv_values, tenor):
	dayCount_maturity = 3 #for the maturity of the swaption
	num_swaptions = len(nv_values) #for a smarter notation

	BPV = np.zeros(num_swaptions); #inizialization
	N_alpha_omega_in_t0 = np.zeros(num_swaptions)
	C_alpha_omega = np.zeros(num_swaptions)

	count_set = range(0,num_swaptions)
	for i in count_set:
		BPV[i] = deltas_fixed_leg[i::] .dot( (discount_factors_semiannual[2*i+2::2] / \
 		discount_factors_semiannual[2*i]) )

		N_alpha_omega_in_t0[i] = 1 - discount_factors_semiannual[len(discount_factors_semiannual)-1] / \
		discount_factors_semiannual[2*i] + \
		(discount_factors_semiannual[2*i:len(discount_factors_semiannual-1)-1] / \
		discount_factors_semiannual[2*i]) .dot ( forward_beta_semiannual[2*i::]-1 )

	fwd_swap_rate = N_alpha_omega_in_t0 / BPV
	dates_semi_aux = dates_semiannual[0:len(dates_semiannual)-2:2]

	dt_t0_alpha = dm.yearfrac([settlement_swaption]*len(dates_semi_aux), dates_semi_aux, 3)
	C_alpha_omega = ( 1 - 1 / (  np.power(1+fwd_swap_rate,tenor)  )  ) / fwd_swap_rate

	prices = discount_factors_semiannual[0:len(discount_factors_semiannual)-2:2] * C_alpha_omega * nv_values * \
	np.sqrt( dt_t0_alpha/(2*np.pi) )

	return [prices, fwd_swap_rate]

#------------------------------------------------------------
#	functions for 'Jamshidian trick'
#------------------------------------------------------------

def f1(x,deltas_fixed_leg,B_alpha_j_fwd,varsigma_alpha_i,i):
	return( deltas_fixed_leg[i::] .dot ( B_alpha_j_fwd * np.exp(-varsigma_alpha_i[1::]*x-0.5*(varsigma_alpha_i[1::]**2)) ) )

def f2(x, B_alphaprime_iota_fwd,varsigma_alphaprime_iota):
	return( B_alphaprime_iota_fwd[1::] .dot (np.exp(-varsigma_alphaprime_iota[1::]*x - 0.5*(varsigma_alphaprime_iota[1::]**2)) ) )

def f3(x, forward_beta_semiannual, B_alphaprime_iota_fwd, nu_alphaprime_iota, i):
	lenB = len(B_alphaprime_iota_fwd)
	return ( forward_beta_semiannual[2*i::] .dot ( B_alphaprime_iota_fwd[0:lenB-1] * np.exp(-nu_alphaprime_iota*x-\
                                                                                              0.5*(nu_alphaprime_iota**2)) ) )
def f_jam(x,deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd, varsigma_alphaprime_iota,
         forward_beta_semiannual, nu_alphaprime_iota, i, strike):

	return( strike*f1(x,deltas_fixed_leg,B_alpha_j_fwd,varsigma_alpha_i,i)+f2(x, B_alphaprime_iota_fwd,varsigma_alphaprime_iota)-\
          f3(x, forward_beta_semiannual, B_alphaprime_iota_fwd, nu_alphaprime_iota, i) )

def S_fun(x,deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd, varsigma_alphaprime_iota,
         forward_beta_semiannual, nu_alphaprime_iota, i):

	return( (f3(x, forward_beta_semiannual, B_alphaprime_iota_fwd, nu_alphaprime_iota, i) -\
             f2(x, B_alphaprime_iota_fwd,varsigma_alphaprime_iota) ) / \
             f1(x,deltas_fixed_leg,B_alpha_j_fwd,varsigma_alpha_i,i) )

def C_alpha_omega(x,deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd, varsigma_alphaprime_iota,
               forward_beta_semiannual, nu_alphaprime_iota, i, tenor):

	return( (1-1./( (1+ S_fun (x,deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd,
                               varsigma_alphaprime_iota, forward_beta_semiannual, nu_alphaprime_iota, i) )**tenor) ) / \
           S_fun(x,deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd, varsigma_alphaprime_iota,
               forward_beta_semiannual, nu_alphaprime_iota, i) )

    
#NB passa solo il tenor giusto
def jam_integrand(x,deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd, varsigma_alphaprime_iota,
         forward_beta_semiannual, nu_alphaprime_iota, i, tenor, strike):

	return( np.exp((-x**2)/2) / np.sqrt(2*np.pi) * \
           max(strike-S_fun(x,deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd, varsigma_alphaprime_iota,
               forward_beta_semiannual, nu_alphaprime_iota, i),0) *\
           C_alpha_omega(x,deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd, varsigma_alphaprime_iota,
               forward_beta_semiannual, nu_alphaprime_iota, i, tenor)
           )



#------------------------------------------------------------
#   Build struct: import data about corporate bonds
#------------------------------------------------------------

def buildStruct():

    # Data for BNPP bonds
    BNPP_dates=list(map(lambda x: datetime.strptime(x, '%d-%b-%Y').toordinal(),['27-Nov-2017','12-Mar-2018','21-Nov-2018','28-Jan-2019',\
    '23-Aug-2019','13-Jan-2021','24-Oct-2022','20-May-2024']))

    BNPP_coupons=list(map(lambda x: x/100,[2.875, 1.5, 1.375, 2,\
     2.5, 2.25, 2.875, 2.375]))

    BNPP_clean=list(map(lambda x: x/100,[105.575, 102.768, 102.555, 104.536, 106.927,\
    106.083,110.281,106.007]));

    BNPP_surv_probs=[1-1.27e-4,1-1.80e-4] #for 2w and 2m 

    BNPP={"dates": BNPP_dates , "coupons":BNPP_coupons,\
    "clean_prices": BNPP_clean, "surv_probs":BNPP_surv_probs}

    # Data for Santander bonds
    Santander_dates=list(map(lambda x: datetime.strptime(x, '%d-%b-%Y').toordinal(),['27-Mar-2017','04-Oct-2017','15-Jan-2018','20-Apr-2018',\
    '14-Jan-2019','13-Jan-2020','24-Jan-2020','14-Jan-2022','10-Mar-2025']))
    
    Santander_coupons=list(map(lambda x:  x/100,[4.000,4.125, 1.750, 0.625, 2.000, 0.875, 4.000,\
    1.125,1.125]))

    Santander_clean=list(map(lambda x: x/100,[105.372,107.358,102.766,\
        99.885,103.984,99.500,112.836,98.166,93.261]))

    Santander_surv_probs=[1-5.42e-4,1-7.73e-4] #for 2w and 2m 

    Santander={"dates": Santander_dates , "coupons":Santander_coupons,\
    "clean_prices": Santander_clean, "surv_probs":Santander_surv_probs}

    return [BNPP,Santander]

#------------------------------------------------------------
#   create vector of coupon payment dates
#------------------------------------------------------------

def create_vector_payment_dates(dates_bonds_maturity, settle, ttl):
	"""Build a unique vector containing all coupon payment dates, and return 
	indexes in order to handle it.
	__________________________________________________________________________
	 INPUT
	 - datesSet:       vector containing all the maturity dates for the bonds;
	 - settle:         settlement date for the basket of bonds;
	 - ttl:            time-to-liquidate.
	--------------------------------------------------------------------------
	 OUTPUT
	 - coupon_payment_dates:  vector that contains all the (unsorted) coupon
	                          payment dates of the bonds (even the ones
	                          between the settlement date and the ttl);
	 - idx1:                  vector that contains the number of coupons 
	                          paid for each bond. 
	                          NB: for practical use, the first element is
	                          set to 0, so the first bond will pay idx(2)
	                          coupons, the second one will pay idx(3) coupons
	                          and so on.
	 - idx1_cs:               cumulative sum of idx_1. This means that
	                          payment dates of the i-th coupon will be
	                          coupon_payment_dates(idx1_cs(i)+1:idx1_cs(i+1))
	                          (indeed if we have N bonds, length(idx1_cs)=N); 
	 - idx2:                  vector that contains the number of coupons paid
	                          between the settlement date and the ttl. The
	                          notation used is the same as idx1 and idx1_cs,
	                          with the first element idx(1) set to 0.
	                          [In our case, since each bond pays annual coupon
	                          and we consider ttl=2weeks or ttl=2monts, idx2
	                          is a binary variable that holds 1 if the coupon 
	                          is paid in between and 0 otherwise].
	--------------------------------------------------------------------------
	 Functions used: dateMoveVec.
	--------------------------------------------------------------------------"""

 # Initial settings
	n_bonds=len(dates_bonds_maturity) # number of bonds
	num_coupon_payments= list(map(lambda x: int(x),np.floor(np.array(dm.yearfrac([settle]*n_bonds,dates_bonds_maturity,0))).tolist()))
	coupon_payment_dates_list= [0]*n_bonds # initialization
	idx1=[0]*(n_bonds+1); #+1 for smarter use in the future
	idx2=[0]*(n_bonds+1); #+1 to use the same notation
 
	# Create vector and indexes
	for k in range(0,n_bonds):
	    # we compute coupon payment dates
	    coupon_payment_dates_list[k]=[0]*(num_coupon_payments[k]+1)

	    for j in range(0,num_coupon_payments[k]+1):
	        coupon_payment_dates_list[k][j]= dm.dateMoveVec_MF([dates_bonds_maturity[k]],\
	                                                        'y', -num_coupon_payments[k]+j, dm.eurCalendar())[0]

	    # we store the number of coupons that will be paid for the i-th bond

	    idx1[k+1]=len(coupon_payment_dates_list[k])

	    # we want to know how many of these coupons fall between settlement and
	    # ttl (in our model, at most one, in the worst cases)
	    idx2[k+1]=len( list(np.extract(np.array(coupon_payment_dates_list[k])<=ttl,\
	        np.array(coupon_payment_dates_list[k]))) )
	#doing the cumulative sum, as explained in the header
	idx1_cs=list(np.cumsum(np.array(idx1))-1)
	
	# converting the dynamic array in a static array
	coupon_payment_dates= [item for sublist in coupon_payment_dates_list\
	                       for item in sublist]
	return [coupon_payment_dates, idx1, idx1_cs, idx2]





#------------------------------------------------------------
#  Compute pi ul
#------------------------------------------------------------
def compute_pi_ul(settle,ttl,payment_dates,idx1,idx1_cs,a,sigma):
	
	# deltas
	dayCount_Act365=3; #day count Act/365 since it a time to maturity

	delta_ttl_payments = dm.yearfrac(ttl*len(payment_dates),payment_dates,dayCount_Act365) #npArray
	delta_settle_ttl = dm.yearfrac([settle],ttl,dayCount_Act365) #npArray

	# cumulated volatility
	ZETA = (sigma/a)*(1-np.exp(-a*delta_ttl_payments)) #npArray
	Sigma_square_tau = (ZETA**2)*( 1-np.exp(-2*a*delta_settle_ttl) ) / (2*a)

	# for each bond, replicating the sigma_square_tau of the maturity
	Sigma_square_N = np.repeat(Sigma_square_tau[idx1_cs[1::]], idx1[1::]) #npArray

	# pi_u
	from scipy.stats import norm

	pi_u = (4+Sigma_square_tau)/2* norm.cdf(np.sqrt(Sigma_square_tau)/2) \
	+ np.sqrt(Sigma_square_tau/(2*np.pi)) * np.exp(-Sigma_square_tau/8)#npArray

	#pi_l

	def f1(eta):
	    r=np.exp(-Sigma_square_N/8)/(np.pi*np.sqrt(1-eta)*np.sqrt(eta))*\
	    np.exp(-eta/2*np.sqrt(Sigma_square_tau)*(np.sqrt(Sigma_square_tau)\
	                                             -np.sqrt(Sigma_square_N) ))
	    return r

	def f2(eta):
	    r=1+np.sqrt(np.pi*(1-eta)/2)*np.sqrt(Sigma_square_N)*\
	    np.exp((1-eta)/8*Sigma_square_N)*norm.cdf(np.sqrt(1-eta)/2*\
	                                              np.sqrt(Sigma_square_N))
	    return r

	def f3(eta):
	    r=1+np.sqrt(np.pi*eta/2)*(2*np.sqrt(Sigma_square_tau)\
	        -np.sqrt(Sigma_square_N))*np.exp(eta/8*\
	        (2*np.sqrt(Sigma_square_tau)-np.sqrt(Sigma_square_N))**2)*\
	        norm.cdf(np.sqrt(eta)/2*(2*np.sqrt(Sigma_square_tau)-np.sqrt(Sigma_square_N)))
	    return r

	def integrand(eta):
	    r=f1(eta)*f2(eta)*f3(eta)
	    return r

	import scipy.integrate as integrate

	pi_l=np.zeros(len(integrand(0.5)))
	for k in range(0,len(pi_l)):
	    pi_l[k]=integrate.quad(lambda eta: integrand(eta)[k],0,1)[0]

	return [list(pi_u),list(pi_l)]

#------------------------------------------------------------
#   create system of equation for Z-bootstrap
#------------------------------------------------------------

def build_system_of_eqn(Z_unknown, flag, coupon_scheme, deltas, B_bar_known, B_interp_known, dirtyPrice, Z_spread_last,
                        delta_last):
	fun = np.zeros(flag)
	fun[0] = -dirtyPrice + np.array(coupon_scheme[0:len(coupon_scheme)-flag]).dot(B_bar_known) + np.array(coupon_scheme[len(coupon_scheme)-flag::])\
	.dot( B_interp_known * np.exp( -deltas * Z_unknown[0:flag] ) )

	fun[1:flag] = (Z_unknown[flag-1]- Z_spread_last) / (deltas[flag-1]-delta_last) * (deltas[0:flag-1]-delta_last) +\
	Z_spread_last - Z_unknown[0:flag-1]
    
	return fun

#------------------------------------------------------------
#	bootstrap_Z_curve
#------------------------------------------------------------

def bootstrap_Z_curve(coupon_payment_dates, idx1, idx1_cs,coupon_on_dates, maturities, discounts_EONIA, dates_EONIA, dirtyPrices):

	settle = dates_EONIA[0]
	dayCount_act365 = 3

	deltas = dm.yearfrac([settle]*len(coupon_payment_dates), coupon_payment_dates, dayCount_act365)

	B_bar = np.zeros(idx1_cs[len(idx1_cs)-1]+1) #initialization
	Z_spread = np.zeros(idx1_cs[len(idx1_cs)-1]+1)

	B_interp = ut.interp_discount(dates_EONIA, discounts_EONIA, coupon_payment_dates)

	def eq_z_cnst(z_cnst):
		return ( -dirtyPrices[0] + np.array(coupon_on_dates[ 0:idx1[1] ]).dot ( B_interp[ 0:idx1[1] ] * 
                                                              np.exp(-z_cnst*deltas[ 0:idx1[1] ]) )  )
	from scipy.optimize import fsolve
	Z_spread[ 0:idx1[1] ] = fsolve(eq_z_cnst,0.005);


	B_bar[0:idx1[1]] = B_interp[0:idx1[1]] *  np.exp(-deltas[0:idx1[1]]*Z_spread[0:idx1[1]]) 

	count_set = range(1,len(maturities))
	for i in count_set:
		fd = idx1_cs[i]+1; #first date (even if before ttl)
		ld = idx1_cs[i+1]; #last date (maturity)
		flag = sum( np.array(coupon_payment_dates[fd:ld+1]) > maturities[i-1] ) #number of unknowns

		idx_sorting = np.argsort( np.insert( deltas[0:fd], 0, 0) )
		deltas_sorted = np.insert( deltas[0:fd], 0, 0)
		deltas_sorted = deltas_sorted[idx_sorting]
		zeta_sorted = np.insert( Z_spread[0:fd], 0, Z_spread[0])
		zeta_sorted = zeta_sorted[idx_sorting]

		Z_spread[fd:ld-flag+1] = np.interp( deltas[fd:ld-flag+1], deltas_sorted, zeta_sorted )
		B_bar[fd:ld-flag+1] = B_interp[fd:ld-flag+1] * ( np.exp(-deltas[fd:ld-flag+1]*Z_spread[fd:ld-flag+1]) )

		if flag == 1:
			B_bar[ld] = (dirtyPrices[i] - np.array(coupon_on_dates[fd:ld-flag+1]).dot(B_bar[fd:ld-flag+1])) / coupon_on_dates[ld]
			Z_spread[ld] = -1/deltas[ld] * np.log( B_bar[ld]/B_interp[ld] )

		else:
			z0 = 1e-2*np.ones(flag)
			Z_spread[ld-flag+1:ld+1] = fsolve(build_system_of_eqn, z0, args=( flag, coupon_on_dates[fd:ld+1], 
                                                              deltas[ld-flag+1:ld+1], B_bar[fd:ld-flag+1], 
                                                              B_interp[ld-flag+1:ld+1], dirtyPrices[i], 
                                                              Z_spread[idx1_cs[i]], deltas[idx1_cs[i]] ) )

			B_bar[ld-flag+1:ld+1] = B_interp[ld-flag+1:ld+1] * np.exp(-deltas[ld-flag+1:ld+1] * Z_spread[ld-flag+1:ld+1])

	return B_bar


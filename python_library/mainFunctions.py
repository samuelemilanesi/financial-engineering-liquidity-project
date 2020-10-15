"""------------------------------------------------------------
# LIBRARY OF MAIN FUNCTIONS
#	- bootstrap_OIS
#	- bootstrap_pseudo
#	- calibrate_model
#	- liquidity_spread_bounds
#	- bond_yield
------------------------------------------------------------"""

#------------------------------------------------------------
#	Used libraries
#------------------------------------------------------------
import pandas as pd
import numpy as np
import scipy as sci
from datetime import datetime

import utility as ut
import dates_manipulation as dm
import auxiliars as aux 

#------------------------------------------------------------
#	Bootstrap EONIA discounts 	
#------------------------------------------------------------

def bootstrap_OIS(dates, rates):
 # Initial settings
	day_count_act360 = 2  
	cleaned_dates_OIS = dates['ois']
	rates_OIS = rates['ois']
	discounts = np.zeros( len(cleaned_dates_OIS) ) #inizialize np vector

 # Dates and Deltas
	deltas_settlement_to_date = dm.yearfrac(dates['settle']*np.ones( len(cleaned_dates_OIS) ), 
		cleaned_dates_OIS, day_count_act360)

	# date 1y after settlement date (int)
	date_1y_after_settlement = dm.dateMoveVec_MF([dates['settle']],'y',1,dm.eurCalendar())

	# index of the first date after 1y from settle
	oneyear = int( sum(( np.ones(len(cleaned_dates_OIS))*\
		[np.array(cleaned_dates_OIS) <= date_1y_after_settlement] ).T) ) 
	
	# delta(t_i,t_i+1)
	deltas = dm.yearfrac( cleaned_dates_OIS[0:len(cleaned_dates_OIS)-1], cleaned_dates_OIS[1::], day_count_act360) 

	# delta(settlement, t_i) for t_i>=1y  NON TROPPO SICURO
	delta_k = np.zeros(len(cleaned_dates_OIS))
	delta_k[oneyear-1] = dm.yearfrac([ dates['settle'] ], [ cleaned_dates_OIS[oneyear-1] ], day_count_act360)
	delta_k[oneyear::] = deltas[oneyear-1::]

	# dates in output
	dates_EONIA = np.zeros(len(cleaned_dates_OIS)+1)
	dates_EONIA[0] = dates['settle']
	dates_EONIA[1::] =  cleaned_dates_OIS

 # Discounts 
	# discounts from OIS up to 1y
	discounts[0:oneyear] = 1 / (1 + deltas_settlement_to_date[0:oneyear]*rates_OIS[0:oneyear]) 
	# bootstrap other discounts 
	for i in range( oneyear-1,len(cleaned_dates_OIS) ):
	    discounts[i] = ( 1 - rates_OIS[i] * delta_k[oneyear-1:i].dot(discounts[oneyear-1:i]) ) / \
	    ( 1+delta_k[i]*rates_OIS[i] )

	# discounts in output
	discounts_EONIA = np.zeros(len(cleaned_dates_OIS)+1)
	discounts_EONIA[0] = 1
	discounts_EONIA[1::] = discounts

	return [dates_EONIA,discounts_EONIA]

#------------------------------------------------------------
#	Bootstrap pseudo-discount curve 	
#------------------------------------------------------------

def bootstrap_pseudo(dates, rates, dates_EONIA, discounts_EONIA):
 # Initial settings
	day_count_act360 = 2
	day_count_30_360 = 6
	day_count_act_365 = 3
 
 # Deltas and dates
 	# delta(settle,6m)
	settle_to_6m_date = dm.yearfrac([ dates['settle'] ],  dates['euribor'] , day_count_act360)

	# delta(6m,1y)
	delta_6m_to_1y = dm.yearfrac( dates['euribor'], [dates['fra_end'][5]] , day_count_act360)

	# delta(fra_start, fra_end) for every fra
	delta_fra = dm.yearfrac(dates['fra_start'][0:5], dates['fra_end'][0:5], day_count_act360);

	# floating dates
	floating_dates = dm.tenordates_MF([dates['settle']], 2*len(dates['swaps']), 'm', 6, dm.eurCalendar());

	# delta(t_i,t_i+1) on fixed dates
	dt_fixed_dates = dm.yearfrac( [dates['settle']] + floating_dates[1:len(floating_dates)-1:2] ,\
                             floating_dates[1:len(floating_dates)+1:2],day_count_30_360 );

 # Discounts up to 1y
	# pseudo-discount at 6m
	discount_pseudo_6m = 1 / (1 + settle_to_6m_date * rates['euribor6m'])

	# pseudo-discount at 1y
	discount_pseudo_1y = discount_pseudo_6m/(1 + delta_6m_to_1y * rates['fra'][5])

	# pseudo-discount 7m to 11m
	discount_7m_to_11m = ut.interp_discount( ut.removenest([ dates['settle'] , dates['euribor'] , dates['fra_end'][5] ]), 
                                         [1,discount_pseudo_6m,discount_pseudo_1y], dates['fra_end'][0:5] )
	# pseudo-discount 1m to 5m
	discount_1m_to_5m = (1 + delta_fra*rates['fra'][0:5]) * discount_7m_to_11m;

	# discount EONIA every semester
	discounts_EONIA_semiannual = ut.interp_discount(dates_EONIA, discounts_EONIA, floating_dates);

	pseudo_discounts_semiannual = np.zeros(len(discounts_EONIA_semiannual))
	pseudo_discounts_semiannual[0]=discount_pseudo_6m
	pseudo_discounts_semiannual[1]=discount_pseudo_1y

 # Define system equation for root_finder
	def equations_1(p):
	    yf_t0_ti = dm.yearfrac([dates['settle']], [floating_dates[2*k+1]], day_count_act_365)
	    yf_t0_timinus1 = dm.yearfrac([dates['settle']], [floating_dates[2*k]], day_count_act_365)
	    yf_t0_timinus2 = dm.yearfrac([dates['settle']], [floating_dates[2*k-1]], day_count_act_365)
	    yf_timinus2_timinus1 = dm.yearfrac([floating_dates[2*k-1]], [floating_dates[2*k]], day_count_act_365)
	    yf_timinus1_ti = dm.yearfrac([floating_dates[2*k]], [floating_dates[2*k+1]], day_count_act_365)
	    yf_timinus2_ti = dm.yearfrac([floating_dates[2*k-1]], [floating_dates[2*k+1]], day_count_act_365)
	    z_timinus2 = -np.log(pseudo_discounts_semiannual[2*k-1])/ yf_t0_timinus2
	    
	    return ( discounts_EONIA_semiannual[2*k] * (pseudo_discounts_semiannual[2*k-1]/p[0]-1) +\
	                discounts_EONIA_semiannual [2*k+1] *(p[0]/p[1]-1) + known_term - I, \
	                (-np.log(p[0])/yf_t0_timinus1 - z_timinus2)/yf_timinus2_timinus1 -\
	                (-np.log(p[1])/yf_t0_ti - z_timinus2)/yf_timinus2_ti)

	# import fsolve in order to solve the system
	from scipy.optimize import fsolve

	for k in range( 1, len(dates['swaps']) ):
	    BPV = discounts_EONIA_semiannual[1:2*k+2:2] .dot(dt_fixed_dates[0:k+1])
	    I = BPV*rates['swaps'][k]
	    
	    forward_euribor6m_semiannual = np.array(  ut.removenest([1, pseudo_discounts_semiannual[0:2*k-1].tolist()])  ) / \
	        np.array( pseudo_discounts_semiannual[0:2*k].tolist() ) -1;
	    
	    known_term = discounts_EONIA_semiannual[0:2*k] .dot(forward_euribor6m_semiannual);
	    
	    x, y =  fsolve(equations_1, (0.5,0.5))
	    pseudo_discounts_semiannual[2*k] = x
	    pseudo_discounts_semiannual[2*k+1] = y

	dates_pseudo = ut.removenest( [ dates['settle'], dates['fra_start'][0:5], dates['euribor'], dates['fra_end'][0:5],\
                dates['fra_end'][5], floating_dates[2::] ] )
	discount_pseudo = [1] + discount_1m_to_5m.tolist() + discount_pseudo_6m.tolist() + discount_7m_to_11m.tolist() + \
                discount_pseudo_1y.tolist() + pseudo_discounts_semiannual[2::].tolist()
	return [dates_pseudo,discount_pseudo]


#------------------------------------------------------------
#	Model Calibration	
#------------------------------------------------------------

def calibrate_model(normal_vols, dates_EONIA, discounts_EONIA, dates_pseudo, discount_pseudo):
 # Initial settings
	tenor = normal_vols['tenor'] # tenors of swaptions
	expiry = normal_vols['expiry']  # swaptions' expiry
	nv_values = normal_vols['values'] 
	settlement = int(dates_EONIA[0]) # settlement date (int)
	dayCount_fixed = 6

 # Deltas and dates
 	# date 1y after the settlement (int)
	date_1y_after_settlement = dm.dateMoveVec_MF([settlement], 'y', 1, dm.eurCalendar())
	# list of semiannual dates 
	dates_semiannual = date_1y_after_settlement+dm.tenordates_MF(date_1y_after_settlement, 2*len(expiry), 'm', 6, dm.eurCalendar())
	
	# delta(ti,t_i+1) on fixed dates
	deltas_fixed_leg = dm.yearfrac(dates_semiannual[0:len(dates_semiannual):2], dates_semiannual[2::2], dayCount_fixed)

 # Discounts 
 	# semiannual discount factors (EONIA and pseudo-discounts)
	discount_factors_semiannual = ut.interp_discount(dates_EONIA, discounts_EONIA, dates_semiannual)
	pseudo_discount_factors_semiannual = ut.interp_discount(dates_pseudo, discount_pseudo, dates_semiannual)
	
	# semiannual forward dicount factors (EONIA and pseudo-discounts) from t_i to t_i+1
	discount_factors_semiannual_fwd = discount_factors_semiannual[1::] / \
    discount_factors_semiannual[0:len(discount_factors_semiannual)-1]
	pseudo_discount_factors_semiannual_fwd = pseudo_discount_factors_semiannual[1::] / \
    pseudo_discount_factors_semiannual[0:len(discount_factors_semiannual)-1]

    # forward beta from t_i to t_i+1
	forward_beta_semiannual = discount_factors_semiannual_fwd / pseudo_discount_factors_semiannual_fwd

 # Compute market prices using Bachelier formula (the function is in auxiliars library)
	[RS_mkt_prices, strike]=aux.normal_vols_swaption_MHJM(dates_semiannual, deltas_fixed_leg,\
	 discount_factors_semiannual, forward_beta_semiannual, settlement, nv_values, tenor)

 #	Define the model-price for reciver cash-settlement swaptions
	def model_prices_swapt(p):
	    #convention: p=[a;sigma;gamma];
	    dayCount_maturity = 3 #Act/365
	    numswapt = len(tenor) #number of swaptions

	    swaption_prices_model  = np.zeros(numswapt) #initialization

	    delta_t0_talpha = dm.yearfrac([settlement]*numswapt, dates_semiannual[0:(2*numswapt-1):2], dayCount_maturity)
	    zeta_alpha_square_vec = p[1]**2 * (1-np.exp(-2*p[0]*delta_t0_talpha)) / (2*p[0])

	    count_set = range(0,numswapt)
	    i=8
	    ld=len(dates_semiannual)

	    for i in count_set:
	        delta_alphaprime_iota = dm.yearfrac([dates_semiannual[2*i]]*(ld-2*i), dates_semiannual[2*i::], 
	                                            dayCount_maturity)
	        v_alphaprime_iota = np.sqrt(zeta_alpha_square_vec[i]) / p[0] *( 1-np.exp(-p[0]*delta_alphaprime_iota) )

	        delta_alpha_i = dm.yearfrac([dates_semiannual[2*i]] * int( (ld+1)/2-i ), dates_semiannual[2*i::2], 
	                                    dayCount_maturity)
	        v_alpha_i = np.sqrt(zeta_alpha_square_vec[i])/p[0]*( 1-np.exp(-p[0] * delta_alpha_i) )
	        nu_alphaprime_iota = v_alphaprime_iota[0:len(v_alphaprime_iota)-1] - p[2] * v_alphaprime_iota[1::]

	        varsigma_alphaprime_iota = (1-p[2]) * v_alphaprime_iota
	        varsigma_alpha_i = (1-p[2]) * v_alpha_i

	        B_alpha_j_fwd = discount_factors_semiannual[2*i+2::2] / discount_factors_semiannual[2*i]
	        B_alphaprime_iota_fwd = discount_factors_semiannual[2*i::] / discount_factors_semiannual[2*i]
	    
	        from scipy.optimize import fsolve

	        x_star = fsolve( aux.f_jam, 0, args=(deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, B_alphaprime_iota_fwd, 
	                                         varsigma_alphaprime_iota,forward_beta_semiannual, nu_alphaprime_iota, i, 
	                                         strike[i]) )
	        x_star=-0.017451279157548
	        from scipy.integrate import quad

	        (x,y) = quad(aux.jam_integrand, -10, x_star, args=(deltas_fixed_leg, B_alpha_j_fwd, varsigma_alpha_i, 
	                                                     B_alphaprime_iota_fwd, varsigma_alphaprime_iota,
	                                                     forward_beta_semiannual, nu_alphaprime_iota, i, tenor[i], 
	                                                     strike[i]) )
	        swaption_prices_model [i] = x * discount_factors_semiannual[2*i]

	    return swaption_prices_model






 # Define the distance function that we have to minimize in order to calibrate the model
	def model_market_distance(p):
	#convention: p=[a;sigma;gamma];
		prices_model = model_prices_swapt(p)
		return( (RS_mkt_prices-prices_model) .dot (RS_mkt_prices-prices_model) )

 # Perform the optimization of model-market distance 
 	# and finding the minimum
	bnds = ( (0, 1), (0, 1), (0,1) ) # defining bounds for parameters

	# The optimal starting point has been found performing a stochastic approach
	x0 = [0.108767495200880,   0.018587703751741,   0.000559789728685] 
	import warnings
	warnings.filterwarnings('ignore', 'The iteration is not making good progress')
	p_opt = sci.optimize.minimize(model_market_distance, x0, bounds=bnds)

 # Calibrated model prices
	RS_calibrated_prices = model_prices_swapt(p_opt.x)

	return[RS_calibrated_prices, RS_mkt_prices, p_opt.x]



#------------------------------------------------------------
#	Liquidity spread bounds	
#------------------------------------------------------------

def liquidity_spread_bounds(bondsInfos, surv_probs, dates_bonds_maturity, calibrated_params,\
 settle, ttl, dirty_prices, discounts_EONIA, dates_EONIA):
	""" Compute upper and lower bounds for the illiquidity price 
	#__________________________________________________________________________
	# INPUT
	# - a,sigma:            calibrated parameters for MHJM model;
	# - ttl:                time to liquidate (date in number);
	# - datesSet:           vector of bond maturities dates;
	# - settle:             settlement date;
	# - coupons:            annual coupons with day-count convention Act/Act
	#                       (one for each bond);
	# - dirtyPrices:        end-of-day mid prices + accrual;
	# - surv_probs:         survival probability for the ttl;
	# - discounts_EONIA:	discount factors for EONIA curve (the first 
	#                       element is 1).
	# - dates_EONIA:        corresponding EONIA dates (the first one is 
	#                       the settlement date);
	#--------------------------------------------------------------------------
	# OUTPUT
	# - pi_u:  vector of upper bound (one for each bond) for illiquidity price;
	# - pi_l:  vector of lower bound (one for each bond) for illiquidity price.
	#--------------------------------------------------------------------------
	# Functions used:   create_vector_payment_dates, compute_pi_ul, 
	#                   bootstrap_Z_curve.
	#__________________________________________________________________________"""

	# variables for more readable code
	#input: a, sigma, ttl, dates_bonds_maturity, settle,\
	#coupons dirty_prices surv_probs discount_EONIA dates_EONIA
	dates_bonds_maturity= bondsInfos["dates"]
	coupons=bondsInfos["coupons"]
	a=calibrated_params[0]
	sigma=calibrated_params[1]
	# Initial settings
	n_bonds=len(dates_bonds_maturity)
	[payment_dates, idx1, idx1_cs, idx2]=aux.create_vector_payment_dates(dates_bonds_maturity, settle, ttl)
	
	coupon_on_dates = list(np.repeat( np.array(coupons), np.array(idx1[1::])))
	
	for k in idx1_cs[1::]:
	    coupon_on_dates[k] = coupon_on_dates[k]+1
	[pi_u, pi_l] = aux.compute_pi_ul(settle, ttl, payment_dates, idx1,\
	                                 idx1_cs, a, sigma)

	#Compute B_bar
	#B_bar = bootstrap_Z_curve(payment_dates, idx1, idx1_cs, ...
	#coupon_on_dates, dates#Set, discount_EONIA, dates_EONIA, dirtyPrices)
	B_bar = aux.bootstrap_Z_curve(payment_dates, idx1, idx1_cs,coupon_on_dates, dates_bonds_maturity, discounts_EONIA, dates_EONIA, dirty_prices)

	# Upper & lower bound

	#for the first time, we do not take into account the coupons that fall
	# between the settle and the ttl, because of 'coupon stripping' (cf. 
	# Baviera-Nassigh-Nastasi, p.10)

	#index of coupons that fall between settlement date and ttl

	tmp=(np.array(idx1_cs[0:len(idx1_cs)-1])+1)*np.array(idx2[1::])
	useless=tmp[list(np.nonzero(tmp)[0])]

	upper_bound =[0]*n_bonds #inizialize
	lower_bound = [0]*n_bonds#inizialize

	addends_u = B_bar*np.array(coupon_on_dates)*(np.array(pi_u)-surv_probs)
	addends_l = B_bar*np.array(coupon_on_dates)*(np.array(pi_l)-surv_probs)

	idx1_cs[0]=-1
	for k in range(0,n_bonds):
	   upper_bound[k] = sum(addends_u[idx1_cs[k]+1:idx1_cs[k+1]+1])
	   lower_bound[k] = sum(addends_l[idx1_cs[k]+1:idx1_cs[k+1]+1])

	return [upper_bound, lower_bound]



#------------------------------------------------------------
#	Bond yield	
#------------------------------------------------------------
def bond_yield(ttl, dates_bonds_maturity, settle,dirty_prices, liquidity_spread, coupons):
	# Creating the vector of coupon payment dates

	n_bonds = len(dates_bonds_maturity)#number of bonds
	[coupon_payment_dates, idx1, idx1_cs, idx2]=aux.create_vector_payment_dates(dates_bonds_maturity, settle, ttl)

	# we consider all the coupons (and no longer the ones after ttl)
	coupon_on_dates = np.repeat( np.array(coupons), np.array(idx1[1::]) )
	coupon_on_dates[ idx1_cs[1::]]= coupon_on_dates[ idx1_cs[1::]] + 1

	# Coupons of interest and deltas
	dayCount_Act365 = 3
	deltas = dm.yearfrac([settle]*len(coupon_on_dates), coupon_payment_dates, dayCount_Act365)

	# Computing P_2bar
	P_2bar = np.array(dirty_prices) - np.array(liquidity_spread)

	# Yield of the liquid bond
	y_T = np.zeros(len(dates_bonds_maturity)) # inizialization
	M = np.zeros(len(dates_bonds_maturity))  # inizialization

	from scipy.optimize import fmin

	for k in range(0,n_bonds):
	    f = lambda y: abs(dirty_prices[k] - np.array(coupon_on_dates[idx1_cs[k]+1:idx1_cs[k+1]+1]).dot(\
	                np.exp(-y*deltas[idx1_cs[k]+1:idx1_cs[k+1]+1])))
	    y_T[k] = fmin(f,0.5,xtol=1e-15,ftol=1e-15,maxiter=100000, maxfun=1e8,disp=False)
	    # we define M to avoid numerical issues
	    g=lambda v: abs(P_2bar[k] - np.array(coupon_on_dates[idx1_cs[k]+1:idx1_cs[k+1]+1]).dot(\
	                    np.exp(-v*deltas[idx1_cs[k]+1:idx1_cs[k+1]+1])))
	    M[k] = fmin(g,0.5,xtol=1e-15,ftol=1e-15,maxiter=100000, maxfun=1e8, disp=False)
	    

	liquidity_yield = M-y_T; #computing liquidity yield

	return [liquidity_yield, y_T]

 
"""------------------------------------------------------------
# LIBRARY OF USEFUL FUNCTIONS:
#	- read_file_excel
#	- yearfrac
#	- removenest
#	- interp_discount
#	- do_plot_error
#	- do_plot_yield
------------------------------------------------------------"""


#------------------------------------------------------------
#	Used libraries
#------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

import dates_manipulation as dm

#------------------------------------------------------------
#	Excel file reader functions	
#------------------------------------------------------------

def read_file_excel(dataFileName):
#Import excel as dataframes
	df_OIS=pd.read_excel(dataFileName, sheet_name="OIS") 
	df_FRA=pd.read_excel(dataFileName, sheet_name="FRA")
	df_vols=pd.read_excel(dataFileName, sheet_name="Vols")
	df_euribor=pd.read_excel(dataFileName, sheet_name="LIBORvs6M")

	#Create dictionaries for quantities of interest

	rates={	"ois": df_OIS[df_OIS.columns[4]].values/100, 
			"fra":df_FRA[df_FRA.columns[6]].values/100,
			"swaps": df_euribor.price.iloc[1:13].values /100,
			"euribor6m": df_euribor.price.iloc[0] /100 }

	dates={	"settle": df_OIS[df_OIS.columns[7]][3].toordinal(),
			"ois": df_OIS[df_OIS.columns[1]].apply(lambda x: to_datetime(x).toordinal()).values,
			"fra_start": df_FRA[df_FRA.columns[2]].apply(lambda x: to_datetime(x).toordinal()).values,
			"fra_end": df_FRA[df_FRA.columns[3]].apply(lambda x: to_datetime(x).toordinal()).values,
			"euribor": [df_euribor[df_euribor.columns[1]][0].toordinal()],
			"swaps": df_euribor[df_euribor.columns[1]].iloc[1:13].apply(lambda x: to_datetime(x).toordinal()).values}

	normal_vols={	"expiry": df_vols[df_vols.columns[1]][1:10].values,
					"tenor": df_vols[df_vols.columns[2]][1:10].values,
					"values": df_vols[df_vols.columns[3]][1:10].values/10000}

	return [rates, dates, normal_vols]

def to_datetime(date):
	timestamp=((date-np.datetime64('1970-01-01T00:00'))/np.timedelta64(1, 's'))
	return datetime.utcfromtimestamp(timestamp)

#------------------------------------------------------------
#	Function removenest
#------------------------------------------------------------
    
def removenest(li):
    return sum(([x] if not isinstance(x, list) else removenest(x)
                for x in li), [])

#------------------------------------------------------------
#	Function interp_discount
#------------------------------------------------------------

def interp_discount(discount_dates,discounts,datesknown):
	day_count_act365=3
	t =  dm.yearfrac([discount_dates[0]]*len(datesknown), datesknown, day_count_act365)  
	tt = dm.yearfrac([discount_dates[0]]*(len(discount_dates)-1), discount_dates[1::], day_count_act365)
	zero_rates = -np.log(discounts[1::]).T/tt 
	count_set = range(0,len(t))
    
	B = np.zeros(len(t))
	tt = tt.tolist()
	zero_rates = removenest( zero_rates.tolist() )
	for i in count_set:
		r = np.interp(t[i], tt, zero_rates)
		B[i] = np.exp(-r*t[i])
	return B

#------------------------------------------------------------
#	dirty_from_clean(dates, coupon, clean_price)
#------------------------------------------------------------

def dirty_from_clean(dates, coupon, clean_prices, settle):
	len_dates=len(dates)
	dayCount_actact=0;
	backward_steps=list(map(int,list(np.ceil(np.array(dm.yearfrac([settle]*len_dates,dates,dayCount_actact))))))
	previous_payment_date = [0]*len_dates
	following_payment_date = [0]*len_dates

	for i in range(0,len_dates):
	    previous_payment_date[i] = dm.dateMoveVec_MF([dates[i]],'y',-backward_steps[i],dm.eurCalendar())[0]
	    following_payment_date[i] = dm.dateMoveVec_MF([dates[i]],'y',-backward_steps[i]+1,dm.eurCalendar())[0]
	ref_yearfrac = dm.yearfrac(previous_payment_date,following_payment_date,0)
	delta = dm.yearfrac(previous_payment_date,[settle]*len_dates,0)

	accrual = np.array(delta)/np.array(ref_yearfrac)*np.array(coupon)
	dirty_prices=list(accrual+np.array(clean_prices))
	return dirty_prices
	

#------------------------------------------------------------
#	do_plot_error
#------------------------------------------------------------

def do_plot_error(settle,dates_bonds_maturity,error2w,error2m,flag):
	#Initial settings
	dayCount_Act365=3 #day count=Act/365 for time to maturity
	maturities=dm.yearfrac([settle]*len(dates_bonds_maturity),\
	                    dates_bonds_maturity,dayCount_Act365)

	# plot

	# title
	if flag=='BNPP':
		plt.title('BNPP price error')
	else:
		plt.title('Santander price error')

	plt.plot(maturities, error2w, 'bd-', maturities, error2m, 'rs:')
	plt.xlabel('Maturity (years)')
	plt.ylabel('Error')
	plt.legend(["ttl=2w","ttl=2m"])
	plt.grid()
	plt.show()

#------------------------------------------------------------
#	do_plot_error
#------------------------------------------------------------

def do_plot_yield(settle,dates_bonds_maturity,bond_yield,liq_yield_2w,liq_yield_2m,flag):
	#Initial settings
	dayCount_Act365=3 #day count=Act/365 for time to maturity
	maturities=dm.yearfrac([settle]*(len(dates_bonds_maturity)),\
		dates_bonds_maturity[0:len(bond_yield)],dayCount_Act365)
	# plot

	# title
	if flag=='BNPP':
		plt.title('BNPP bond yields')
	else:
		plt.title('Santander bond yields')

	plt.plot( np.array(maturities), 10000*np.array(bond_yield), 'bs-',\
	         np.array(maturities), 10000*(np.array(bond_yield)+liq_yield_2w), 'rs:',\
	        np.array(maturities), 10000*(np.array(bond_yield)+liq_yield_2m), 'gs:') 
	plt.xlabel('Maturity (years)')
	plt.ylabel('yields (bps)')
	plt.legend(["liquid","ttl = 2w","ttl = 2m"])
	plt.grid()
	plt.show()
"""------------------------------------------------------------
# LIBRARY OF FUNCTIONS ABOUT DATES MANIPULATION:
#	- eurCalendar
#	- yearfrac
#	- isbusday
#	- busday_next
#	- busday_prev
#	- first_date_MF
#	- eom
#	- dateadd
#	- dateMoveVec_MF
#	- tenordates
------------------------------------------------------------"""


#------------------------------------------------------------
#	Used libraries
#------------------------------------------------------------
import pandas as pd
import numpy as np
from datetime import datetime

import utility as ut

#------------------------------------------------------------
# 	eurCalendar
#------------------------------------------------------------

def eurCalendar():
    h=[730120, 730231, 730234, 730241, 730479, 730480, 730486, 730588, 730591, 730606, 730844, 
            730844, 730845, 730851, 730938, 730941, 730971, 731209, 731210, 731216, 731323, 731326, 
            731326, 731336, 731574, 731575, 731581, 731680, 731683, 731702, 731940, 731941, 731947, 
            731947, 732030, 732033, 732067, 732305, 732306, 732312, 732415, 732418, 732432, 732670, 
            732670, 732671, 732677, 732772, 732775, 732797, 733035, 733036, 733042, 733122, 733125, 
            733125, 733163, 733401, 733402, 733408, 733507, 733510, 733528, 733766, 733767, 733773, 
            733773, 733864, 733867, 733893, 734131, 734132, 734138, 734249, 734252, 734258, 734496, 
            734496, 734497, 734503, 734599, 734602, 734624, 734862, 734863, 734869, 734956, 734959, 
            734959, 734989, 735227, 735228, 735234, 735341, 735344, 735354, 735592, 735593, 735599, 
            735599, 735691, 735694, 735719, 735957, 735958, 735964, 736048, 736051, 736085, 736323, 
            736323, 736324, 736330, 736433, 736436, 736450, 736688, 736689, 736695, 736783, 736786, 
            736786, 736815, 737053, 737054, 737060, 737168, 737171, 737180, 737418, 737419, 737425, 
            737425, 737525, 737528, 737546, 737784, 737785, 737791, 737882, 737885, 737911, 738149, 
            738149, 738150, 738156, 738260, 738263, 738276, 738514, 738515, 738521, 738617, 738620, 
            738620, 738641, 738879, 738880, 738886, 738974, 738977, 739007, 739245, 739246, 739252,
            739252, 739359, 739362, 739372, 739610, 739611, 739617, 739709, 739712, 739737, 739975, 
            739975, 739976, 739982, 740066, 740069, 740102, 740340, 740341, 740347, 740451, 740454, 
            740454, 740468, 740706, 740707, 740713, 740801, 740804, 740833, 741071, 741072, 741078, 
            741078, 741186, 741189, 741198, 741436, 741437, 741443, 741543, 741546, 741563, 741801, 
            741801, 741802, 741808, 741893, 741896, 741929, 742167, 742168, 742174, 742278, 742281, 
            742281, 742294, 742532, 742533, 742539, 742635, 742638, 742659, 742897, 742898, 742904, 
            742904, 742985, 742988, 743024, 743262, 743263, 743269, 743370, 743373, 743390, 743628,
            743628, 743629, 743635, 743727, 743730, 743755, 743993, 743994, 744000, 744112, 744115, 
            744115, 744120, 744358, 744359, 744365, 744462, 744465, 744485, 744723, 744724, 744730, 
            744730, 744819, 744822, 744851, 745089, 745090, 745096, 745204, 745207, 745216, 745454, 
            745454, 745455, 745461, 745554, 745557, 745581, 745819, 745820, 745826, 745911, 745914, 
            745914, 745946, 746184, 746185, 746191, 746296, 746299, 746312, 746550, 746551, 746557, 
            746557, 746653, 746656, 746677, 746915, 746916, 746922, 747003, 747006, 747042, 747280, 
            747280, 747281, 747287, 747388, 747391, 747407, 747645, 747646, 747652, 747745, 747748, 
            747748, 747773, 748011, 748012, 748018, 748123, 748126, 748138, 748376, 748377, 748383, 
            748383, 748480, 748483, 748503, 748741, 748742, 748748, 748837, 748840, 748868, 749106, 
            749106, 749107, 749113, 749222, 749225, 749234, 749472, 749473, 749479, 749572, 749575, 
            749575, 749599, 749837, 749838, 749844, 749929, 749932, 749964, 750202, 750203, 750209, 
            750209, 750314, 750317, 750329, 750567, 750568, 750574, 750664, 750667, 750695, 750933, 
            750933, 750934, 750940, 751049, 751052, 751060, 751298, 751299, 751305, 751406, 751409, 
            751409, 751425, 751663, 751664, 751670, 751756, 751759, 751790, 752028, 752029, 752035, 
            752035, 752141, 752144, 752156, 752394, 752395, 752401, 752498, 752501, 752521, 752759, 
            752759, 752760, 752766, 752848, 752851, 752886, 753124, 753125, 753131, 753233, 753236, 
            753236, 753251, 753489, 753490, 753496, 753590, 753593, 753617, 753855, 753856, 753862, 
            753862, 753947, 753950, 753982, 754220, 754221
            ]
    return h

#------------------------------------------------------------
#	Function yearfrac
#------------------------------------------------------------
#dates1 and dates2 must have the same length
def yearfrac(dates1, dates2, idx):
	len_dates2 = len(dates2)
	count1 = range(0,len_dates2)
	yf=np.zeros(len(dates2))
	for i in count1:
#ACT/ACT
		if idx==0:
			from datetime import datetime
			dt1 = datetime.fromordinal(dates1[i])
			import datetime
			next_year = datetime.date(dt1.year+1, dt1.month, dt1.day).toordinal()
			yf[i] = (dates2[i]-dates1[i])/(next_year-dates1[i])
#ACT/360
		if idx==2:
			yf[i] = (dates2[i]-dates1[i])/360
#ACT/365
		if idx==3:
			yf[i] = (dates2[i]-dates1[i])/365
#30/360
		if idx==6:
			from datetime import datetime
			dt1 = datetime.fromordinal(dates1[i])
			dt2 = datetime.fromordinal(dates2[i])
			if dt1.day == 31:
				dt1.day == 30
			if dt2.day == 31:
				dt2.day == 30
			yf[i] = (360*(dt2.year-dt1.year)+30*(dt2.month-dt1.month)+(dt2.day-dt1.day)) / 360

	return yf

#------------------------------------------------------------
#	Function isbusday
#------------------------------------------------------------

def isbusday(date_1, calendar):
	from datetime import datetime
	date_2 = datetime.fromordinal(date_1)
	fest = calendar
	if date_2.weekday() > 4: #saturdays and sundays
		return 0
	if date_1 in fest:
		return 0
	return 1

#------------------------------------------------------------
#	Function busday_next
#find the first business day after date_1 (which is not date_1!)
#------------------------------------------------------------

def busday_next(date_1, calendar):
	flag = 0  
	while flag == 0:
		date_1 += 1
		flag = isbusday(date_1, calendar)
	return date_1

#------------------------------------------------------------
#	Function busday_prev
#find the first business day before date_1 (which is not date_1!)
#------------------------------------------------------------

def busday_prev(date_1, calendar):
	flag = 0  
	while flag == 0:
		date_1 -= 1
		flag = isbusday(date_1, calendar)
	return date_1   

#------------------------------------------------------------
#	Function first_date_MF
#------------------------------------------------------------

def first_date_MF(date_vec, calendar):
	len_vec = len(date_vec)
	count_set = range(0,len_vec)
	dates_ok = [0]*len_vec #return an integer
	for i in count_set:
		if isbusday(date_vec[i], calendar) == 0:
			date_tmp = busday_next(date_vec[i], calendar)
			from datetime import datetime
			date_ext_1 = datetime.fromordinal(date_vec[i])
			from datetime import datetime
			date_ext_2 = datetime.fromordinal(date_tmp)
			if date_ext_1.month != date_ext_2.month:    
				dates_ok[i] = busday_prev(date_vec[i], calendar)
			else:
				dates_ok[i] = date_tmp
		else:
			dates_ok[i] = date_vec[i]
	return  dates_ok

#------------------------------------------------------------
#	Function eom (end of month)
#------------------------------------------------------------

def eom(year, month):
    thirty = [4, 6, 9, 11]
    thirtyone = [1, 3, 5, 7, 8, 10, 12]
    if month in thirtyone:
        return 31
    if month in thirty:
        return 30
    if year%4 == 0:
        return 29
    return 28

#------------------------------------------------------------
#	Function dateadd
#------------------------------------------------------------

def dateadd(startdate, datepart, num): #scalar input
	if datepart == 'd':
		fdate = startdate + num
		return fdate
	if datepart == 'w':
		fdate = startdate + 7*num
		return fdate
	if datepart == 'm':
		from datetime import datetime
		dt1 = datetime.fromordinal(startdate)
		dt2 = np.zeros(3)
		dt2[0] = dt1.year + np.floor( (dt1.month+num-1)/12 )
		dt2[1] = (dt1.month + num -1)%12 +1
		dt2[2] = min(dt1.day, eom(dt2[0],dt2[1]))
		import datetime
		fdate = datetime.date( int(dt2[0]), int(dt2[1]), int(dt2[2]) ).toordinal()
		return fdate
	if datepart == 'y':
		from datetime import datetime
		dt1 = datetime.fromordinal(startdate)
		import datetime
		fdate = datetime.date( dt1.year+num, dt1.month, dt1.day ).toordinal()
		return fdate
    
#------------------------------------------------------------
#	Function dateMoveVec_MF
#------------------------------------------------------------

def dateMoveVec_MF(startdates, datepart, num, calendar): #passa sempre una lista come input startdates!
	len_startdates = len(startdates)
	if len_startdates==1:
		temp_dates_moved=dateadd(startdates[0] , datepart, num)
		dates_moved = first_date_MF([temp_dates_moved],calendar)
		return dates_moved
	count_set = range( 0 , len_startdates )
	dates_moved = [0] * len_startdates
	for i in count_set:
		dates_moved[i] = dateadd(startdates[i], datepart, num)
	dates_moved = first_date_MF(dates_moved, calendar)
	return dates_moved

#------------------------------------------------------------
#	Function tenordates
#------------------------------------------------------------

def tenordates_MF(startdates, nsteps, datepart, num, calendar):
	Dates = [0]*nsteps
	count_set = range(0,nsteps)
	for i in count_set:
		Dates[i] = dateMoveVec_MF(startdates, datepart, num*(i+1), calendar)
	return ut.removenest(Dates)
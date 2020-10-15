% CORPORATE BOND LIQUIDITY: exericises 3 and 4
% Financial Engineering: Politecnico Milano
%
% In order to run the script:
% >> runFinalProject_Group6_3_4
%
% Last Modified: 14.06.2019
%
% REMARK: we assume to have a face value equal to 1, so everything is
%         rescaled according to this assumption

clc
clear all
close all
addpath(genpath('.'));

%% Discount bootstrap (as in point i)
formatData='dd/mm/yyyy'; % modified in order to use the same format date of the spreadsheet 
[datesSet_dirty, ratesSet, normal_vols] = readExcelData('curves20150910_project.xlsx',formatData);

datesSet=clean_date(datesSet_dirty); %we want to use business dates only

% Perform bootstrap for EONIA curve
[dates_EONIA,discounts_EONIA] = bootstrap_OIS(datesSet, ratesSet); 

%% iii) Upper and lower bounds of liquidity spread
% Given data
[datesSet_2, coupon, clean_price, surv_probs] = buildStruct();
a = 0.1294; sigma = 0.0126; % Calibrated parameters of point ii)

% time to liquidate: 2 cases
two_weeks = dateMoveVec(datesSet.settlement, 'w', 2, 'MF', eurCalendar);
two_months = dateMoveVec(datesSet.settlement, 'm', 2, 'MF', eurCalendar);

%% Dirty prices from clean prices (we add the accrual term)

dirty_prices_BNPP = dirty_from_clean(datesSet_2.BNPP, ...
    coupon.BNPP, clean_price.BNPP, datesSet.settlement);
dirty_prices_Santander = dirty_from_clean(datesSet_2.Santander, ...
    coupon.Santander, clean_price.Santander, datesSet.settlement);
 
%% BNPP
% 2weeks              
[upper_bound_BNPP_2w, lower_bound_BNPP_2w]=liquidity_spread_bounds...
    (a, sigma, two_weeks, datesSet_2.BNPP, datesSet.settlement, ...
     coupon.BNPP, dirty_prices_BNPP, surv_probs.BNPP(1), ...
     discounts_EONIA, dates_EONIA);
error_BNPP_2w=abs(upper_bound_BNPP_2w-lower_bound_BNPP_2w);

% 2months            
[upper_bound_BNPP_2m, lower_bound_BNPP_2m]=liquidity_spread_bounds...
    (a, sigma, two_months, datesSet_2.BNPP, datesSet.settlement,...
     coupon.BNPP, dirty_prices_BNPP, surv_probs.BNPP(2),...
     discounts_EONIA, dates_EONIA);
error_BNPP_2m=abs(upper_bound_BNPP_2m-lower_bound_BNPP_2m);

do_plot_error(datesSet.settlement, datesSet_2.BNPP, error_BNPP_2w, ...
    error_BNPP_2m,1) % plot BNPP

%% Santander 
% 2weeks

[upper_bound_Santander_2w, lower_bound_Santander_2w]=liquidity_spread_bounds...
    (a, sigma, two_weeks, datesSet_2.Santander, datesSet.settlement, ...
     coupon.Santander, dirty_prices_Santander, surv_probs.Santander(1), ...
     discounts_EONIA, dates_EONIA);
error_Santander_2w=abs(upper_bound_Santander_2w-lower_bound_Santander_2w);

% 2months

[upper_bound_Santander_2m, lower_bound_Santander_2m]=liquidity_spread_bounds...
    (a, sigma, two_months, datesSet_2.Santander, datesSet.settlement,...
     coupon.Santander, dirty_prices_Santander, surv_probs.Santander(2),...
     discounts_EONIA, dates_EONIA);
error_Santander_2m=abs(upper_bound_Santander_2m-lower_bound_Santander_2m);

do_plot_error(datesSet.settlement, datesSet_2.Santander, error_Santander_2w,...
    error_Santander_2m, 2) % plot Santander
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Point iv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BNPP

% 2 weeks
[liq_yield_BNPP_2w, bond_yield_BNPP_2w]=bond_yield...
    (two_weeks,datesSet_2.BNPP, datesSet.settlement,...
    dirty_prices_BNPP, upper_bound_BNPP_2w, coupon.BNPP);
% 2 months
[liq_yield_BNPP_2m, bond_yield_BNPP_2m]=bond_yield...
    (two_months,datesSet_2.BNPP, datesSet.settlement,...
    dirty_prices_BNPP, upper_bound_BNPP_2m, coupon.BNPP);

do_plot_yield(datesSet.settlement, datesSet_2.BNPP, bond_yield_BNPP_2w,...
    liq_yield_BNPP_2w,liq_yield_BNPP_2m, 1)% plot BNPP

%% Santander
% 2 weeks
[liq_yield_Santander_2w,bond_yield_Santander_2w]=bond_yield...
    (two_weeks,datesSet_2.Santander, datesSet.settlement,...
    dirty_prices_Santander, upper_bound_Santander_2w, coupon.Santander);

% 2 months
[liq_yield_Santander_2m,bond_yield_Santander_2m]=bond_yield...
    (two_months,datesSet_2.Santander, datesSet.settlement,...
    dirty_prices_Santander, upper_bound_Santander_2m, coupon.Santander);

do_plot_yield(datesSet.settlement, datesSet_2.Santander, bond_yield_Santander_2w,...
    liq_yield_Santander_2w,liq_yield_Santander_2m, 2) % plot Santander

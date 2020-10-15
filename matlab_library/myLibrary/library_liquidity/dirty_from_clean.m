function dirty_price=dirty_from_clean(datesSet, coupon, clean_price, settle)

% Compute dirty prices from clean prices.
%__________________________________________________________________________
% INPUT
% - datesSet:      vector containing the maturities of the bonds;
% - coupon:        vector of coupons (one for each bond);
% - clean_price:   vector of clean prices of the bonds;
% - settle:        settlement date.
%--------------------------------------------------------------------------
% OUTPUT
% - dirty_price:    vector of corresponding dirty prices (clean prices +
%                   accrual).
%--------------------------------------------------------------------------
% Functions used: dateMoveVec.
%--------------------------------------------------------------------------
% Last Modified: 04.06.2019
%__________________________________________________________________________

%% Dates of interest
% first we find how many years we have to go backwards to find the first
% coupon payment before the settlement (for each bond)
backward_steps=ceil(yearfrac(settle,datesSet));

% then we find the date of the first coupon payment before the settlement
% (for each bond)
previous_payment_date = zeros(length(datesSet),1);
following_payment_date = zeros(length(datesSet),1);

 for i=1:length(datesSet)
    previous_payment_date(i) = ...
       dateMoveVec(datesSet(i),'y',-backward_steps(i), 'MF',eurCalendar);
   	following_payment_date(i) = ...
       dateMoveVec(datesSet(i),'y',-backward_steps(i)+1,'MF',eurCalendar);
 end

%% Accrual: we compute the accrual according to the dates computed before 
% computing time interval between one payment date and the following one
ref_yearfrac = yearfrac(previous_payment_date,following_payment_date);

%time interval for which we have to 'give back' a part of the coupon
delta = yearfrac(previous_payment_date,settle); 

% Accrual to be paid in addition to the clean price
accrual = delta ./ ref_yearfrac .* coupon; 
% accrual=delta .* coupon;

%% Dirty price

dirty_price=accrual+clean_price;

end
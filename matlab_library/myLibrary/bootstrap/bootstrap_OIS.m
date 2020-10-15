function [dates_EONIA,discounts_EONIA] = bootstrap_OIS(datesSet, ratesSet)

% Perform bootstrap for the EONIA curve, given a dataset of Overnight 
% Indexed Swaps (OIS).
%__________________________________________________________________________
% INPUT
% - datesSet:       struct containing the following elements:
%                       - settlement:   settlement date;
%                       - OIS:          vector of expiry dates of the 
%                                       Overnight Indexed Swaps;
% - ratesSet:       struct containing the following elements:
%                       - OIS:          vector of the corresponding
%                                       quoted OIS rates.
%--------------------------------------------------------------------------
% OUTPUT
% - dates_EONIA:        reset dates for OIS (the same as input datesSet,
%                       but the first one is the settlement date);
% - discounts_EONIA:    discount factors for EONIA curve in the 
%                       corresponding dates (the first element is 1).
%--------------------------------------------------------------------------
% Last Modified: 06.06.2019
%__________________________________________________________________________


%     
% REMARK: we suppose that after the first year, we have 
% yearly datas for the OIS
%

%% Initial settings

day_count_act360 = 2;

dates = datesSet.OIS; %in order to have a more readable code
discounts = zeros(length(dates),1); %inizialization

% deltas from t_0 to t_i
deltas_settlement_to_date = ...
    yearfrac(datesSet.settlement, dates, day_count_act360);

%% If OIS lasts less than one year

% finding the 1y OIS (then we need to use the other formula)
date_1y_after_settlement = ...
    dateMoveVec(datesSet.settlement,'y',1,'MF',eurCalendar);
%index of the last maturity that does not exceed 1y
oneyear = sum(dates <= date_1y_after_settlement); 

% building the discount curve
discounts(1:oneyear) = 1 ./ ...
    (1 + deltas_settlement_to_date(1:oneyear).*ratesSet.OIS(1:oneyear));

%% If the OIS has maturity longer than one year

% deltas from t_i to t_{i+1}
deltas = yearfrac( datesSet.OIS(1:end-1), datesSet.OIS(2:end), ...
             day_count_act360);

% constructing a vector of deltas
delta_k = [deltas(1:oneyear-1);
           yearfrac(datesSet.settlement, datesSet.OIS(oneyear), ...
                day_count_act360);
           deltas(oneyear:end)];
%
% REMARK: in practice, delta_k under the first year are not used!
% This is just to have a smarter notation for the indexes.
%

% completing the discount curve
for i = oneyear+1:length(dates)
    discounts(i) = ( 1 - ratesSet.OIS(i) * delta_k(oneyear:i-1)' * ...
        discounts(oneyear:i-1) ) / ( 1+delta_k(i).*ratesSet.OIS(i) );
end

%% Output

%including settlement date
dates_EONIA     = [datesSet.settlement; dates];
discounts_EONIA = [1; discounts];

end
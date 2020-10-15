function  B  = interp_discounts( discount_dates, discounts, dates)

% Interpolate discount factors using a linear-on-zero-rates interpolation.
% 'interp_discounts' has the same input/output structure as 'interp1'.
%__________________________________________________________________________
% INPUT:
% - discount_dates:     dates in which the discount curve is known;
% - discounts:          corresponding discount factors;
% - dates:              dates in which we want to compute discount factors.
%--------------------------------------------------------------------------
% OUTPUT:
% - B:                  discount factors in 'dates'.
%--------------------------------------------------------------------------
% Last Modified: 07.06.2019
%__________________________________________________________________________


%% Settings 
day_count_act365=3; % zero rates day-count convention

%% Interpolating

% dates in which we need the discount factors
t = yearfrac(discount_dates(1), dates, day_count_act365);

% dates in which we know the discount factors
tt = yearfrac(discount_dates(1), discount_dates(2:end), day_count_act365); 

% known zero rates
zero_rates = -log(discounts(2:end))./tt;

% linar interpolation on zero rates 
r = interp1(tt, zero_rates, t);

%% Output

% from zero rates to discount factors
B = exp(-r.*t);

end
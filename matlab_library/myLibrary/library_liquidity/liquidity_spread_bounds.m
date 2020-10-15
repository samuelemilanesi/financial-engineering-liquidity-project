function [upper_bound, lower_bound] = liquidity_spread_bounds(a, sigma, ttl,...
          datesSet, settle, coupons, dirtyPrices, surv_probs, discount_EONIA,...
          dates_EONIA)

% Compute upper and lower bounds for the illiquidity price 
%__________________________________________________________________________
% INPUT
% - a,sigma:            calibrated parameters for MHJM model;
% - ttl:                time to liquidate (date in number);
% - datesSet:           vector of bond maturities dates;
% - settle:             settlement date;
% - coupons:            annual coupons with day-count convention Act/Act
%                       (one for each bond);
% - dirtyPrices:        end-of-day mid prices + accrual;
% - surv_probs:         survival probability for the ttl;
% - discounts_EONIA:	discount factors for EONIA curve (the first 
%                       element is 1).
% - dates_EONIA:        corresponding EONIA dates (the first one is 
%                       the settlement date);
%--------------------------------------------------------------------------
% OUTPUT
% - pi_u:  vector of upper bound (one for each bond) for illiquidity price;
% - pi_l:  vector of lower bound (one for each bond) for illiquidity price.
%--------------------------------------------------------------------------
% Functions used:   create_vector_payment_dates, compute_pi_ul, 
%                   bootstrap_Z_curve.
%--------------------------------------------------------------------------
% Last Modified: 09.06.2019
%__________________________________________________________________________

%% Build the vector of dates

n_bonds=length(datesSet);

[payment_dates, idx1, idx1_cs, idx2] = ...
    create_vector_payment_dates(datesSet,settle,ttl);

coupon_on_dates = repelem( coupons, idx1(2:end) );
coupon_on_dates(idx1_cs(2:end)) = coupon_on_dates(idx1_cs(2:end)) + 1;

%% Compute Pi_u and Pi_l

[pi_u, pi_l] = compute_pi_ul(settle, ttl, payment_dates, idx1, ...
    idx1_cs, a, sigma);

%% Compute B_bar

B_bar = bootstrap_Z_curve(payment_dates, idx1, idx1_cs, ...
    coupon_on_dates, datesSet, discount_EONIA, dates_EONIA, dirtyPrices);

%% Upper & lower bound

% for the first time, we do not take into account the coupons that fall
% between the settle and the ttl, because of 'coupon stripping' (cf. 
% Baviera-Nassigh-Nastasi, p.10) -> we exploit the vector idx2

upper_bound = zeros(n_bonds,1);
lower_bound = zeros(n_bonds,1);

addends_u = B_bar .* coupon_on_dates .* (pi_u-surv_probs) ;
addends_l = B_bar .* coupon_on_dates .* (pi_l-surv_probs) ;

for k=1:n_bonds
   upper_bound(k) = sum( addends_u(idx1_cs(k)+idx2(k+1)+1:idx1_cs(k+1)) );
   lower_bound(k) = sum( addends_l(idx1_cs(k)+idx2(k+1)+1:idx1_cs(k+1)) );
end

end
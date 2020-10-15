function [liquidity_yield, y_T] = bond_yield(ttl, datesSet, settle, ...
    dirtyPrices, liquidity_spread, coupons)

% Compute bond yields for the liquid bond
%__________________________________________________________________________
% INPUT
% - ttl:                time to liquidate;
% - datesSet:           vector of maturities of the bonds;
% - settle:             settlement date;
% - dirtyPrices:        end-of-day mid prices (dirty prices);
% - liquidity_spread:   vector of liquidity_spreads (one for each bond);
% - coupon:             vector of annual coupons with day-count convention 
%                       Act/Act (one for each bond).
%--------------------------------------------------------------------------
% OUTPUT
% - liquidity_yield:    vector of liquidity yields (one for each bond);
% - y_T:                vector of bond yields (one for each bond).
%--------------------------------------------------------------------------
% Functions used: create_vector_payment_dates.
%--------------------------------------------------------------------------
% Last Modified: 08.06.2019
%__________________________________________________________________________

%% Creating the vector of coupon payment dates

n_bonds = length(datesSet); %number of bonds

[payment_dates, idx1, idx1_cs, ~] = ...
    create_vector_payment_dates(datesSet,settle,ttl);

% we consider all the coupons (and no longer the ones after ttl)
coupon_on_dates = repelem( coupons, idx1(2:end) );
coupon_on_dates(idx1_cs(2:end)) = coupon_on_dates(idx1_cs(2:end)) + 1;

%% Coupons of interest and deltas

dayCount_Act365 = 3;

deltas = yearfrac(settle, payment_dates, dayCount_Act365);

%% Computing P_2bar

P_2bar = dirtyPrices - liquidity_spread;

%% Yield of the liquid bond

y_T = zeros(length(datesSet),1); % inizialization
M = zeros(length(datesSet),1);   % inizialization

options = optimset('Display','off','TolFun',1e-15); % for fsolve

for i=1:n_bonds
    f=@(y) dirtyPrices(i) - coupon_on_dates(idx1_cs(i)+1:idx1_cs(i+1))' * ...
            exp( -y*deltas(idx1_cs(i)+1:idx1_cs(i+1)) );
    y_T(i) = fsolve(f, 0.5, options);
    
    %we define M to avoid numerical issues
    g=@(M) P_2bar(i) - coupon_on_dates(idx1_cs(i)+1:idx1_cs(i+1))' * ...
            exp( -M*deltas(idx1_cs(i)+1:idx1_cs(i+1)) );
    M(i) = fsolve(g, 0.5, options);
end

liquidity_yield = M-y_T; %computing liquidity yield

end
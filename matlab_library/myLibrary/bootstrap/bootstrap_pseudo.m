function [dates_pseudo,discount_pseudo] = ...
    bootstrap_pseudo(datesSet,ratesSet, dates_EONIA, discounts_EONIA)

% Perform bootstrap for the EONIA curve, given a dataset of Euribor-based
% financial instruments as described below. The function uses FRAs in order
% to deduce the discounts from 1m to 5m using 'Mr.Crab's algorithm'
% (cf. 'A Note on Dual-Curve Construction: Mr. Crab's Bootstrap', 
% Baviera-Cassaro, 2015)
%__________________________________________________________________________
% INPUT
% - datesSet:       struct containing the following elements:
%                   	- settlement:   settlement date;
%                       - euribor:      expiry of a 6m swap;
%                   	- fra:          matrix of settlement dates and 
%                                       expiries for Forward Rate 
%                                       Agreements (at least 1x7, 2x8, 3x9,
%                                       4x10, 5x11, 6x12);
%                       - swaps:        vector of maturities of the swaps;
%
% - ratesSet:       struct containing the following elements:
%                       - Euribor6m:    swap rate of a 6m swap, i.e. the
%                                       value of the 6m-euribor;
%                   	- fra:          vector of forward rates for FRAs;
%                       - swaps:        vector of swap rates.
%
% - dates_EONIA:        reset dates for OIS (the same as input datesSet,
%                       but the first one is the settlement date);
% - discounts_EONIA:    discount factors for EONIA curve in the 
%                       corresponding dates (the first element is 1).
%--------------------------------------------------------------------------
% OUTPUT
% - dates_pseudo:       reset dates for 6m Euribor pseudo discount curve 
%                       (basically the same as input datesSet, but the 
%                       first one is the settlement date);
% - discount_pseudo:    discount factors for 6m Euribor pseudo discount 
%                       curve corresponding dates (the first element is 1).
%--------------------------------------------------------------------------
% Functions used: interp_discounts, tenordates.
%--------------------------------------------------------------------------
% Last Modified: 06.06.2019
%__________________________________________________________________________

%% Settings
day_count_act360 = 2; % for floating leg 
day_count_30_360 = 6; % for fixed leg 
%% 6 Months pseudo discount
% we use 6m Euribor swap to deduce the 6m pseudo discount factor
settle_to_6m_date = yearfrac(datesSet.settlement, datesSet.euribor, ...
    day_count_act360);
% computing the 6m pseudo discount factor
discount_pseudo_6m = 1 / (1 + settle_to_6m_date * ratesSet.Euribor6m);

%% 1y pseudo discount
% we use 6x12 FRA to compute the 1y pseudo discount factor

delta_6m_to_1y = yearfrac(datesSet.euribor, datesSet.fra(6,2), ...
    day_count_act360);

% computing the 1y pseudo discount factor, exploiting fwd discount factor
discount_pseudo_1y = discount_pseudo_6m/(1+delta_6m_to_1y*ratesSet.fra(6));

%% Mr.Crab
% we use some FRAs to obtain the discounts in the initial buckets

% linear interpolation on zero rate, using 6m and 1y values
discount_7m_to_11m = interp_discounts( ...
        [datesSet.settlement; datesSet.euribor; datesSet.fra(6,2)], ...
        [1; discount_pseudo_6m; discount_pseudo_1y], ...
        datesSet.fra(1:5,2) );

delta_fra = yearfrac(datesSet.fra(1:5,1), datesSet.fra(1:5,2), ...
    day_count_act360);

% computing pseudo discounts via Mr.Crab algorithm
discount_1m_to_5m = (1+delta_fra.*ratesSet.fra(1:5)) .* discount_7m_to_11m;

%% 2y to 12y pseudo discounts
% we use swaps-vs-Euribor6m to compute the 2y to 12y pseudo discounts

% payment dates of the floating leg (every 6m)
floating_dates = tenordates(datesSet.settlement, ...
    2*length(datesSet.swaps), 'm', 6, 'MF', eurCalendar);

% interpolating (linear on zero rates) the EONIA discounts for floating leg
discounts_EONIA_semiannual = interp_discounts(dates_EONIA, ...
    discounts_EONIA, floating_dates);

% year-fraction for fixed leg
dt_fixed_dates = yearfrac( ...
    [datesSet.settlement;floating_dates(2:2:end-2)],...
    floating_dates(2:2:end),day_count_30_360 );

pseudo_discounts_semiannual = zeros(24,1); % inizialization: 24 half-years
pseudo_discounts_semiannual(1:2) = [ discount_pseudo_6m;
                                     discount_pseudo_1y]; %known discounts

for k=2:length(ratesSet.swaps)
    BPV = discounts_EONIA_semiannual(2:2:2*k)' * dt_fixed_dates(1:k);
    I = BPV*ratesSet.swaps(k); %NPV of the fixed leg
    
    % NB: not exactly the forward_euribor6m. Indeed this vector contains
    % the forward_euribor6m multiplied by the corresponding year-fraction
    forward_euribor6m_semiannual = ...
        [1; pseudo_discounts_semiannual(1:2*k-3)] ./ ...
        pseudo_discounts_semiannual(1:2*k-2) -1;

    % known part of the NPV of the fixed leg
    known_term = discounts_EONIA_semiannual(1:2*k-2)' * ...
        forward_euribor6m_semiannual;
   
    %interpolation and roots-finding ('two_eqns' is defined in the footer)
    missing_pseudo_discounts = @(P) two_eqns(P, known_term, I, ...
        pseudo_discounts_semiannual, discounts_EONIA_semiannual, ...
        floating_dates, datesSet.settlement, k);
    x0 = [1;1]; %initial point for search

    options = optimset('Display','off','TolFun',1e-18);
    pseudo_discounts_semiannual(2*k-1:2*k) = ...
        fsolve(missing_pseudo_discounts,x0,options); %solving the system

end

%% Final output

dates_pseudo = [ datesSet.settlement;
                 datesSet.fra(1:5,1);
                 datesSet.euribor;
                 datesSet.fra(1:5,2)
                 datesSet.fra(6,2);
                 floating_dates(3:end)];
        
discount_pseudo = [ 1;
                    discount_1m_to_5m;
                    discount_pseudo_6m;
                    discount_7m_to_11m;
                    discount_pseudo_1y;
                    pseudo_discounts_semiannual(3:end)];

end


%--------------------------------------------------------------------------
% Non-linear system to compute pseudo discounts from swap rate
%--------------------------------------------------------------------------
function fun = two_eqns(P, known_term, I, pseudo_discounts_semiannual,...
                discounts_EONIA_semiannual, floating_dates, settlement, k)

    % imposing NPV=0
    fun(1) = discounts_EONIA_semiannual(2*k-1) * ...
        (pseudo_discounts_semiannual(2*(k-1))/P(1)-1) +...
        discounts_EONIA_semiannual(2*k)*(P(1)./P(2)-1) + known_term - I;
    
    % the half-year discount should fulfill a linear-on-zero-rate
    % interpolation condition
	fun(2) = interp_discounts( ...
             [settlement; floating_dates(1:2*k-2); floating_dates(2*k)], ...
             [1; pseudo_discounts_semiannual(1:2*k-2) ;P(2)], ...
             floating_dates(2*k-1) ...
             ) - P(1);
end
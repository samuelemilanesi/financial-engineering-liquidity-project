function [a, sigma, gamma] = calibrate_model(dates_EONIA, ...
    discounts_EONIA, dates_pseudo, discounts_pseudo, normal_vols)

% Calibrate volatility parameters [a, sigma, gamma] for the Multicurve
% HJM model (MHJM) on ATM Cash-Settlement (CS) diagonal receiver swaptions.
%__________________________________________________________________________
% INPUT
% - dates_EONIA:        vector of dates for EONIA discount curve
%                       (first element should be the settlement date);
% - discounts_EONIA:	vector of corresponding EONIA discount factors;
% - dates_pseudo:       vector of dates for Euribor6m discount curve
%                       (first element should be the settlement date);
% - discounts_pseudo:	vector of corresponding Euribor6m discount factors;
% - normal_vols:        struct containing the following fields (vectors):
%                           - values:   values of normal volatilities;
%                           - expiry:   expiries of the swaptions (i.e.
%                                       settlement dates of the underlying
%                                       swaps);
%                           - tenor:    tenor of the underlying swaps.
%--------------------------------------------------------------------------
% OUTPUT
% - [a, sigma, gamma]:  calibrated parameters for MHJM model.
%--------------------------------------------------------------------------
% Functions used: tenordates, interp_discounts, normal_vols_swaption_MHJM,
%                 compute_model_swaption_MHJM
%--------------------------------------------------------------------------
% Last Modified: 07.06.2019
%__________________________________________________________________________

%% Settings
% renaming the variable in order to have a more readable code
tenor = normal_vols.tenor;
expiry = normal_vols.expiry;
nv_values = normal_vols.values;
settlement = dates_EONIA(1); %dates_EONIA(1) contains the settlement date

% daycounts
dayCount_fixed = 6; %30/360, fixed leg

%% Dates

% finding the 1y date (expiry of the first swaption)
date_1y_after_settlement = dateMoveVec(settlement,'y',1,'MF',eurCalendar);

% dates every 6m from 1y(included) to 10y(included)
dates_semiannual = [  date_1y_after_settlement;
                      tenordates(date_1y_after_settlement, ...
                            2*length(expiry),'m',6,'MF',eurCalendar)  ];

% deltas for the fixed leg (for BPV computation)
deltas_fixed_leg = yearfrac(dates_semiannual(1:2:end-2),...
    dates_semiannual(3:2:end),dayCount_fixed);
                        
%% Market prices from volatilities, part#1: discounts

% deducing discount and pseudo-discount factors for the dates in
% dates_semiannual via linear-on-zero-rates interpolation
discount_factors_semiannual = interp_discounts(dates_EONIA, ...
    discounts_EONIA, dates_semiannual);
pseudo_discount_factors_semiannual = interp_discounts(dates_pseudo, ...
    discounts_pseudo, dates_semiannual);

% computing the 'multiplicative spread', i.e. beta(t_0;t_i,t_{i+1}):
% - step#1: forward semiannual discount factors, i.e. B(t_0;t_i,t_{i+1}),
% where 'i' starts from 1y included (first swaption expiry)
discount_factors_semiannual_fwd = discount_factors_semiannual(2:end)./ ...
    discount_factors_semiannual(1:end-1);
% - step#2: forward semiannual pseudo discount factors, i.e. 
% B_tilde(t_0;t_i,t_{i+1})
pseudo_discount_factors_semiannual_fwd = ...
    pseudo_discount_factors_semiannual(2:end) ./ ...
    pseudo_discount_factors_semiannual(1:end-1);
% - step#3: multiplicative spread, i.e. beta(t_0;t_i,t_{i+1})
forward_beta_semiannual = discount_factors_semiannual_fwd./...
    pseudo_discount_factors_semiannual_fwd;

%% Market prices from volatilities, part#2: strikes and prices

% strikes = forward swap rates (ATM swaptions)
[RS_mkt_prices, strike] = normal_vols_swaption_MHJM( ...
    dates_semiannual, deltas_fixed_leg, discount_factors_semiannual, ...
    forward_beta_semiannual, settlement, nv_values, tenor);

%% Minimisation
% notation: from now on, x=[a;sigma;gamma];

% we define a function handle of 'x' in order to perform the minimisation
model_price_swaption = @(a,sigma,gamma) ...
    compute_model_swaption_MHJM(a, sigma, gamma, settlement, ...
    tenor, strike, dates_semiannual, deltas_fixed_leg, ...
    discount_factors_semiannual, forward_beta_semiannual);

% objective function
d=@(x) norm( model_price_swaption(x(1),x(2),x(3)) - RS_mkt_prices ).^2;

% starting point
x0=[0.108767495200880   0.018587703751741   0.000559789728685]; %best x0

%--------------------------------------------------------------------------
% with fminsearch
%--------------------------------------------------------------------------

t = fminsearch(d,x0);

%--------------------------------------------------------------------------
% % with a global minimum algorithm
%--------------------------------------------------------------------------

% % domain of the parameters
% lb=[0;0;0];
% ub=[0.3;0.1;0.001];
% 
% % Imposing a finer search
% options = optimoptions('fmincon','Algorithm','interior-point');
% options = optimoptions(options, 'StepTolerance', 1e-18, ...
%           'ConstraintTolerance', 1e-18);
%
% problem = createOptimProblem('fmincon', ...
%     'objective',@(x) d(x), ...
%     'x0', x0, ...
%     'lb', lb, 'ub', ub,'options',options);
%
% warning('off','all')
%    
% ms = MultiStart('Display','iter');
% t = run(ms,problem, 50);

%% Defining the output

a       = t(1);
sigma   = t(2);
gamma   = t(3);

%% Plot

figure(3)
plot(RS_mkt_prices,'rs-')
grid on
hold on
plot(model_price_swaption(a,sigma,gamma),'bd:')
xlabel('Expiry')
ylabel('Price of the swaptions')
title('Calibration of the model')
legend('Mkt prices','MHW prices','Location','NorthEast')
ylim([0 0.035])

end
function [pi_u, pi_l] = compute_pi_ul(settle, ttl, payment_dates, idx1, ...
	idx1_cs, a, sigma)

% Compute the vectors of parameters pi_u, pi_l that are used in the
% formulas for upper and lower bounds of illiquidity price (cf. 'A closed 
% formula for illiquid corporate bonds and an application to the European 
% market', Baviera-Nassigh-Nastasi, 2019)
%__________________________________________________________________________
% INPUT
% - settle:             settlement date for the set of bonds;
% - ttl:                time-to-liquidate;
% - payment_dates:      vector that contains all the (unsorted) coupon
%                       payment dates of the bonds (even the ones between
%                       the settlement date and the ttl);
% - idx1, idx1_cs:      indexes for mananaging payment_dates (cf. 
%                       create_vector_payment_dates for furhter info);
% - a, sigma:           calibrated parameters for MHJM model.
%--------------------------------------------------------------------------
% OUTPUT
% - pi_u:  vector of upper parameters (one for each bond);
% - pi_l:  vector of lower parameters (one for each bond).
%--------------------------------------------------------------------------
% Last Modified: 09.06.2019
%__________________________________________________________________________

%% deltas

dayCount_Act365=3; % day count Act/365 since it a time to maturity

delta_ttl_payments = yearfrac(ttl,payment_dates,dayCount_Act365);
delta_settle_ttl = yearfrac(settle,ttl,dayCount_Act365);

%% Cumulated volatility

ZETA = (sigma/a) * ( 1-exp(-a*delta_ttl_payments) );    

Sigma_square_tau = (ZETA.^2) * ( 1-exp(-2*a*delta_settle_ttl) ) / (2*a);

% for each bond, replicating the sigma_square_tau of the maturity
Sigma_square_N = repelem( Sigma_square_tau(idx1_cs(2:end)), idx1(2:end) );

%% pi_u 

pi_u = (4+Sigma_square_tau)/2 .* normcdf(sqrt(Sigma_square_tau)/2) + ...
    sqrt(Sigma_square_tau/(2*pi)) .* exp(-Sigma_square_tau/8);

%% pi_l

f1 = @(eta) exp(-Sigma_square_N/8)./(pi*sqrt(1-eta).*sqrt(eta)).*...
     exp(-eta/2.*sqrt(Sigma_square_tau) .* ( sqrt(Sigma_square_tau) ...
     -sqrt(Sigma_square_N) ));
f2 = @(eta) 1 + sqrt(pi*(1-eta)/2) .* sqrt(Sigma_square_N) .* ...
     exp( (1-eta)/8.*Sigma_square_N ) .* normcdf(sqrt(1-eta)/2 .* ...
     sqrt(Sigma_square_N));
f3 = @(eta) 1 + sqrt(pi*eta/2) * ...
    (2*sqrt(Sigma_square_tau)-sqrt(Sigma_square_N)) .* ...
    exp(eta/8*(2*sqrt(Sigma_square_tau)-sqrt(Sigma_square_N)).^2) .* ...
    normcdf(sqrt(eta)/2*(2*sqrt(Sigma_square_tau)-sqrt(Sigma_square_N)));

integrand = @(eta) f1(eta).*f2(eta).*f3(eta);

pi_l = integral(integrand, 0, 1, 'ArrayValued', true);

end
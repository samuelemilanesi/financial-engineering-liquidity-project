function [B_bar] = bootstrap_Z_curve( payment_dates, idx1, idx1_cs,...
    coupon_on_dates, maturities, discount_EONIA, dates_EONIA, dirtyPrices)

% Perform the bootstrap on Z-spread: Zeta-spread curve is assumed to be
% constant up to the maturity of the bond with the lowest maturity and
% it is linearly interpolated afterwards.
%__________________________________________________________________________
% INPUT
% - payment_dates:          vector that contains all the (unsorted) coupon
%                           payment dates of the bonds (even the ones
%                           between the settlement date and the ttl);
% - idx1, idx1_cs:          indexes for mananaging payment_dates (cf. 
%                           create_vector_payment_dates for furhter info);
% - maturities:             liquid bond maturities;
% - discounts_EONIA:        discount factors for EONIA curve (the first 
%                           element is 1).
% - dates_EONIA:            corresponding EONIA dates (the first one is 
%                           the settlement date);
% - dirtyPrices:            dirty prices of the liquid bonds.
%--------------------------------------------------------------------------
% OUTPUT
% - B_bar:                  vector of defaultable ZC rates, having the
%                           same length as payment_dates.
%--------------------------------------------------------------------------
% Functions used:   interp_discounts, build_system_of_eqn
%--------------------------------------------------------------------------
% Last Modified: 08.06.2019
%__________________________________________________________________________

%% Settings

options = optimset('Display','off','TolFun',1e-12);  %for fsolve

settle = dates_EONIA(1); %settlement date

dayCount_act365 = 3; %daycount convention for z_spread

deltas = yearfrac(settle, payment_dates, dayCount_act365);

B_bar = zeros(idx1_cs(end),1); %inizialization
Z_spread = zeros(idx1_cs(end),1); %inizialization

% interpolating the discount factors (we already have the entire curve)
B_interp = interp_discounts(dates_EONIA, discount_EONIA, payment_dates);

%% First z_spread (constant case)

%solving numerically the first equation (Z(1) is constant)
f=@(Z_start) -dirtyPrices(1) + coupon_on_dates( 1:idx1(2) )' * ...
                ( B_interp(1:idx1(2)) .* exp(-Z_start*deltas(1:idx1(2))) );

Z_spread(1:idx1(2)) = fsolve(f,0.005,options);

% computing the associated B_bars
B_bar(1:idx1(2)) = B_interp(1:idx1(2)) .* ...
                        exp( -deltas(1:idx1(2)) .* Z_spread(1:idx1(2)) );
%% Interpolation of Z-spread curve

% preliminar step: due to interp1 issues, we select the dates that appears
% more than once in payment_dates
[~,filter1] = unique(deltas,'stable');

for i=2:length(maturities)
    fd = idx1_cs(i)+1; %first date (even if before ttl)
    ld = idx1_cs(i+1); %last date (maturity)
    flag = sum( payment_dates(fd:ld) > maturities(i-1) ); %number of unknown
    % escamotage: we put to negative numbers to the deltas for which we 
    % haven't computed the z_spread yet
	deltas_tmp = [deltas(1:fd-1); -(fd:idx1_cs(end))'];
	deltas_interp = deltas_tmp(filter1); %deltas for interpolation
	clear deltas_tmp
    
    
	% interpolating the known z_spreads (the ones that fall before the
	% last available data, i.e. the maturity of the {i-1}-th bond)
    Z_spread(fd:ld-flag) = interp1( [0;deltas_interp], ...
        [Z_spread(1);Z_spread(filter1)], deltas(fd:ld-flag) );
	% deducing the corresponding B_bar
	B_bar(fd:ld-flag) = B_interp(fd:ld-flag) .* ...
        exp( -deltas(fd:ld-flag) .* Z_spread(fd:ld-flag) );
    
	if flag == 1 %caso base (analytic equation)
        
        % deducing the unknown B_bar
        B_bar(ld) = (dirtyPrices(i) - coupon_on_dates(fd:ld-flag)' *...
            B_bar(fd:ld-flag)) / coupon_on_dates(ld);
        
        % and computing the corresponding z_spread (for interpolation)
        Z_spread(ld) = -1/deltas(ld) * log( B_bar(ld)/B_interp(ld) );
	else % at least two unknowns, non-trivial case (numeric solution)
        
        % creating the system of equations that will be solved
        missing_Z_spreads = @(Z_unknown) build_system_of_eqn( Z_unknown, ...
            flag, coupon_on_dates(fd:ld), deltas(ld-flag+1:ld), ...
            B_bar(fd:ld-flag), B_interp(ld-flag+1:ld), dirtyPrices(i), ...
            Z_spread(idx1_cs(i)), deltas(idx1_cs(i)) );
        
        % initial point
        z0 = 1e-2*ones(flag,1);
        
        % solving the system
        Z_spread(ld-flag+1:ld) = fsolve(missing_Z_spreads, z0, options);
        
        % and clearing the variables
        clear missing_Z_spreads;
        clear Z_unknown;
        
        % lastly, compute again the B_bars using the new z-spreads 
        B_bar(ld-flag+1:ld) = B_interp(ld-flag+1:ld) .* ...
            exp(-deltas(ld-flag+1:ld) .* Z_spread(ld-flag+1:ld));

	end
end

end
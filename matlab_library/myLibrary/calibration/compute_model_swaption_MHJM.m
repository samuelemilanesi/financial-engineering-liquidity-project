function prices = compute_model_swaption_MHJM(a, sigma, gamma, ...
    swaptions_settlement, tenor, strike, dates_semiannual, ...
    deltas_fixed_leg, discount_factors_semiannual, forward_beta_semiannual)

% Computes MHJM model price for a complete set (i.e. one per year, same 
% settlement day, first swaption expires at 1y) of ATM Cash-Settlement 
% diagonal reciver swaptions vs Euribor6m (fixed leg paid annually, 
% floating semiannually).
%__________________________________________________________________________
% INPUT
% - a, sigma, gamma:                parameters  for MHJM model;
% - swaptions_settlement:           settlement date for all the swaptions, 
%                                   i.e. t0;
% - tenor:                          vector of tenors of the swaptions;
% - strike:                         vector of strikes of the swaptions,
%                                   i.e. corresponding fwd swap rate, since
%                                   they are all ATM;
% - dates_semiannual:               vector of all the payment dates of the 
%                                   floating legs of the underlying swaps,
%                                   (it should also contain as first 
%                                   element the expiry date of the
%                                   first swaption);
% - deltas_fixed_leg:               vector of year-fractions for the fixed 
%                                   leg of the underlying swaps (it should 
%                                   contain yearly fractions, i.e. 
%                                   [delta(1y,2y); delta(2y,3y); ...
%                                   delta(3y,4y); etc.] with an appropriate
%                                   day-count);
% - discount_factors_semiannual:	vector of EONIA discount factors. It 
%                                   should contain discount factor wrt t0,
%                                   i.e. B(t0;t0,t_i) with 'i' that goes
%                                   from the settlement of the first
%                                   swaption (1y) to the common expiry of
%                                   the underlying swaps (extrema included)
%                                   and with a semiannual time-step;
% - forward_beta_semiannual:        vector of multiplicative spreads
%                                   between EONIA and Euribor6m curve. It 
%                                   should include semiannual spreads
%                                   i.e. beta(t0;t_i,t_{i+1}) with 'i' that 
%                                   runs from the settlement of the first
%                                   swaption (1y) to the date before the
%                                   common expiry of the underlying swaps.
%--------------------------------------------------------------------------
% OUTPUT
% - prices:                         vector of prices of the swaptions,
%                                   accordingly to MHJM model.
%--------------------------------------------------------------------------
% Last Modified: 07.06.2019
%__________________________________________________________________________

%% Initial settings
dayCount_maturity=3; %Act/365
prices=zeros(length(tenor),1); %inizialization

%% Computing the prices

% zeta_alpha_square_vec can be computed before the cycle because it
% depends only on tenor

delta_t0_talpha = yearfrac(swaptions_settlement, ...
    dates_semiannual(1:2:(2*length(tenor)-1)), dayCount_maturity);
zeta_alpha_square_vec = sigma^2 * (1-exp(-2*a*delta_t0_talpha)) / (2*a);

for i = 1:length(tenor)
    % defining the whole set of model parameters for the i-th swaption
    % (2*i-1) always means alpha(i), i.e. the expiry of the i-th swaption
    delta_alphaprime_iota = yearfrac(dates_semiannual(2*i-1), ...
        dates_semiannual(2*i-1:end), dayCount_maturity);
    v_alphaprime_iota = sqrt(zeta_alpha_square_vec(i))/a*( 1-exp(-a* ...
        delta_alphaprime_iota) );
    
    delta_alpha_i = yearfrac(dates_semiannual(2*i-1),...
        dates_semiannual(2*i-1:2:end),dayCount_maturity);
    v_alpha_i = sqrt(zeta_alpha_square_vec(i))/a*( 1-exp(-a*...
        delta_alpha_i) );
    
    nu_alphaprime_iota = v_alphaprime_iota(1:end-1) - ...
        gamma * v_alphaprime_iota(2:end); % one element shorter
    
    varsigma_alphaprime_iota = (1-gamma) * v_alphaprime_iota;
    varsigma_alpha_i = (1-gamma) * v_alpha_i;

    % useful discount factors for first summation
    B_alpha_j_fwd = discount_factors_semiannual(2*i+1:2:end) / ...
        discount_factors_semiannual(2*i-1);
    % useful discount factors for second and third summations
    B_alphaprime_iota_fwd = discount_factors_semiannual(2*i-1:end) / ...
        discount_factors_semiannual(2*i-1);
    
    % functions for 'Jamshidian trick'
    f1=@(x) sum( deltas_fixed_leg(i:end) .* B_alpha_j_fwd .* ...
        exp(-varsigma_alpha_i(2:end)*x-0.5*varsigma_alpha_i(2:end).^2) );
    f2=@(x) sum( B_alphaprime_iota_fwd(2:end) .* ...
            exp(-varsigma_alphaprime_iota(2:end)*x - ...
                0.5*varsigma_alphaprime_iota(2:end).^2) );
    f3=@(x) sum( forward_beta_semiannual(2*i-1:end) .* ...
        B_alphaprime_iota_fwd(1:end-1) .* ...
        exp(-nu_alphaprime_iota*x-0.5*nu_alphaprime_iota.^2) );
    
    % summing up
    f=@(x) strike(i)*f1(x)+f2(x)-f3(x);
    

    
    %finding the ' Jamshidian's x_star '
    options = optimset('Display','off');
    x_star = fsolve(f,0,options);
    
    % obtaining the integrand function
    S = @(x) (f3(x)-f2(x))./f1(x);
    C_alpha_omega=@(x) ( 1-1./((1+S(x)).^tenor(i)) ) ./ S(x);
    integrand = @(x) normpdf(x) .* C_alpha_omega(x) .* max(strike(i)-S(x),0);
    
    %computing the integral
    prices(i) = discount_factors_semiannual(2*i-1) * ...
        quadgk(integrand, -10, x_star);

end

end
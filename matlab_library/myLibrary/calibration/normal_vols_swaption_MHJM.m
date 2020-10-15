function [prices, fwd_swap_rate] = normal_vols_swaption_MHJM( ...
    dates_semiannual, deltas_fixed_leg, discount_factors_semiannual, ...
    forward_beta_semiannual, settlement_swaption, nv_values, tenor)

% Computes market price for a complete set (i.e. one per year, same 
% settlement day, first swaption expires at 1y) of ATM Cash-Settlement 
% diagonal reciver swaptions vs Euribor6m (fixed leg paid annually, 
% floating semiannually), given a corresponding set of normal volatilities. 
% The function uses the normal-Black formula (a.k.a. Bachelier formula) to 
% convert normal volatilities into prices and works under the hypothesis of 
% multicurve HJM framework.
%__________________________________________________________________________
% INPUT
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
%                                   common expiry of the underlying swaps;
% - settlement_swaption:            settlement of the swaptions, i.e. t0;
% - nv_values:                      vector of normal volatilities;
% - tenor:                          vector of tenors of the swaptions.
%--------------------------------------------------------------------------
% OUTPUT
% - prices:                         vector of prices of the swaptions;    
% - fwd_swap_rate:                  forward swap rates of the underlying
%                                   swaps (i.e. the strikes of the 
%                                   swaptions).
%--------------------------------------------------------------------------
% Last Modified: 07.06.2019
%__________________________________________________________________________

%% Settings
dayCount_maturity = 3; %for the maturity of the swaption
num_swaptions = length(nv_values); %for a smarter notation

%% Computing forward swap rates

BPV = zeros(num_swaptions,1); %inizialization
N_alpha_omega_in_t0 = zeros(num_swaptions,1);

for i=1:num_swaptions
    % computing 'fwd' BPV
    BPV(i) = deltas_fixed_leg(i:end)' * ...
     (discount_factors_semiannual(2*i+1:2:end) / ...
      discount_factors_semiannual(2*i-1));
    
    % computing 'fwd' NPV of the floating leg (multicurve case)
    N_alpha_omega_in_t0(i) = 1 - discount_factors_semiannual(end) / ...
        discount_factors_semiannual(2*i-1) + ...
        discount_factors_semiannual(2*i-1:end-1)' ./ ...
        discount_factors_semiannual(2*i-1) * ...
        ( forward_beta_semiannual(2*i-1:end)-1 );
end

% obtaining the forward swap rates, i.e. the strikes (ATM swaptions)
fwd_swap_rate = N_alpha_omega_in_t0./BPV;

%% Bachelier Formula

dt_t0_alpha = yearfrac(settlement_swaption, ...
    dates_semiannual(1:2:end-2), dayCount_maturity);

% equivalent of the BPV for the CS swaptions
C_alpha_omega = ( 1-1./((1+fwd_swap_rate).^tenor) ) ./ fwd_swap_rate;

%Bachelier Formula for CS receiver swaptions
prices = discount_factors_semiannual(1:2:end-2) .* ...
    C_alpha_omega .* nv_values .* sqrt(dt_t0_alpha/(2*pi));

end
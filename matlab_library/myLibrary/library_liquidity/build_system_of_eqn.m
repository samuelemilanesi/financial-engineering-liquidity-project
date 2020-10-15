function fun = build_system_of_eqn(Z_unknown, flag, coupon_scheme, ...
    deltas, B_bar_known, B_interp_known, dirtyPrice, Z_spread_last, ...
    delta_last)

% Auxiliary function that builds the system of equations that are solved
% in the bootstrap_Z_curve. The system is meant to be a non-linear system
% related to a unique bond, the i-th bond.
%__________________________________________________________________________
% INPUT
% - Z_unknown:       unknown z_spreads (in the Z-bootstrap, this is the 
%                    vector of unknown that is used as variable in the 
%                    function  handle);
% - flag:            number of unknown;
% - coupon_scheme:   vector of coupon payments for the i-th bond;
% - deltas:          yearfraction (settle-payment) for the unknown: its
%                    size should be (flag,1);
% - B_bar_known:     known values of B_bar, i.e. for all the payment dates
%                    for which we already know the z-spread: its size
%                    should be ( length(coupon_scheme)-flag,1 );
% - B_interp_known:  EONIA discount factors in the last 'flag' dates, i.e.
%                    the ones for which we do not know the z_spread;
% - dirtyPrice:      dirty price of the i-th bond;
% - delta_last:      largest delta for which we know the z-spread (should
%                    be the maturity of the {i-1}-th bond;
% - Z_spread_last:   corresponding z-spread.
%--------------------------------------------------------------------------
% OUTPUT
% - fun:    vector containing the non-linear system that has to be solved.
%--------------------------------------------------------------------------
% Last Modified: 09.06.2019
%__________________________________________________________________________

%Z_spread last and delta_last are the last known deltas/spread, i.e. the
%ones of the previous maturity
fun = zeros(flag,1);
fun(1) = -dirtyPrice + coupon_scheme(1:end-flag)'*B_bar_known + ...
    coupon_scheme(end-flag+1:end)' * ( B_interp_known .* ...
    exp( -deltas .* Z_unknown(1:flag) ) );

% linear interpolation condition: basically, for each equation we impose 
% that the slopes of the two segments must be the same
fun(2:flag) = (Z_unknown(flag)- Z_spread_last) ./ (deltas(flag)-delta_last) .* ...
    (deltas(1:flag-1)-delta_last) + Z_spread_last - Z_unknown(1:flag-1);
end
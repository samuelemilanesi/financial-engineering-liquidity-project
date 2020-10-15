function []=do_plot_yield(settle, datesSet, bond_yield, ...
                          liq_yield_2w,liq_yield_2m, flag) 
% Plot obligor yields and yields obtained for illiquid bonds with ttl equal
% to two weeks and two months
%__________________________________________________________________________
% INPUT
% - settle:       	 	settlement date (number)
% - datesSet:			vector of bond maturities
% - bond_yield:         vector of bond yields
% - liq_yield_2w:       liquidity yield for 2w ttl (number)
% - liq_yield_2w:       liquidity yield for 2m ttl (number)
% - flag:               1 for BNPP, 2 for Santander
%--------------------------------------------------------------------------
% OUTPUT
% 	no output
%--------------------------------------------------------------------------
% Last Modified: 07.06.2019
%__________________________________________________________________________
%% 
dayCountAct365=3; % day count=Act/365 for time to maturity
maturities=yearfrac(settle,datesSet,dayCountAct365);
%% Plot
figure
plot(maturities, 1e4*bond_yield, 'b-s') %liquid
hold on
plot(maturities, 1e4*(bond_yield+liq_yield_2w), 'r--^')
grid on
plot(maturities, 1e4*(bond_yield+liq_yield_2m), 'g--o')

%% Legend and labels
xlabel('Maturity (years)')
ylabel('yields (bps)')
legend('liquid','\tau =2 weeks','\tau =2 months','location','northwest')

%% title
if flag==1
    title('BNPP bond yields')
else
    title('Santander bond yields')
end

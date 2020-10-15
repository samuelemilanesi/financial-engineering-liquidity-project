function []=do_plot_error(settle, datesSet, error2w, error2m,flag)
% Plot the difference between the upper and lower bounds for the
% illiquidity price Delta_tau for obligor bonds
%__________________________________________________________________________
% INPUT
% - settle:       	 	settlement date (number)
% - datesSet:			vector of bond maturities
% - error2w:            spread upper-lower bound in case 2w (number)
% - error2m:            spread upper-lower bound in case 2m (number)
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
plot(maturities, error2w, 'b-s')
hold on
plot(maturities, error2m, 'r--^')
grid on

%% Legend and labels
xlabel('Maturity (years)')
ylabel('Error')
legend('\tau =2 weeks','\tau =2 months','location','northwest')

%% title
if flag==1
    title('BNPP price error')
else
    title('Santander price error')
end
end
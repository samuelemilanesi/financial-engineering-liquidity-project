function  plot_discount(dates_EONIA, discounts_EONIA, dates_pseudo, ...
                discounts_pseudo,formatData) 
% Plot EONIA and Euribor6m curves for the Multicurve HJM model (MHJM) and
% the corresponding zero rates curves.
%__________________________________________________________________________
% INPUT
% - dates_EONIA:        vector of dates for EONIA discount curve
%                       (first element should be the settlement date);
% - discounts_EONIA:	vector of corresponding EONIA discount factors;
% - dates_pseudo:       vector of dates for Euribor6m discount curve
%                       (first element should be the settlement date);
% - discounts_pseudo:	vector of corresponding Euribor6m discount factors;
% - formatData:         format for the dates that appear in the plot.
%--------------------------------------------------------------------------
% OUTPUT
% 	no output
%--------------------------------------------------------------------------
% Functions used: tenordates, interp_discounts, normal_vols_swaption_MHJM,
%                 compute_model_swaption_MHJM
%--------------------------------------------------------------------------
% Last Modified: 07.06.2019
%__________________________________________________________________________

%% Settings

day_count_maturity=3; %Act/365, for zero rates

%% Computing the zero rates

yf1 = yearfrac(dates_EONIA(1),dates_EONIA(2:end),day_count_maturity);
zero_rates_EONIA = -log(discounts_EONIA(2:end)) ./ yf1;

yf2 = yearfrac(dates_pseudo(1),dates_pseudo(2:end),day_count_maturity);
zero_rates_pseudo = -log(discounts_pseudo(2:end)) ./ yf2;

%% Plot of discount factors

figure(1)
plot(dates_EONIA, discounts_EONIA, 'r-')
grid on
hold on
plot(dates_pseudo, discounts_pseudo, 'b-')

xt=datenum({'14/Sep/2015','14/Sep/2017','14/Sep/2019','14/Sep/2021',...
    '14/Sep/2023','14/Sep/2025','14/Sep/2027'});
xticks(datenum(xt))
xlim([xt(1),xt(end)])
xticklabels(datestr(xt,formatData))
xtickangle(45)
title('Comparison between the two discount curves')
legend('OIS curve','Euribor6m curve','Location','NorthEast')

%% Plot of the corresponding zero rates

figure(2)
plot(dates_EONIA(2:end), zero_rates_EONIA, 'r-')
grid on
hold on
plot(dates_pseudo(2:end), zero_rates_pseudo, 'b-')

xticks(datenum(xt))
xlim([xt(1),xt(end)])
xticklabels(datestr(xt,formatData))
xtickangle(45)
title('Comparison between the zero rates')
legend('OIS zero rates','Euribor6m zero rates','Location','NorthEast')

end
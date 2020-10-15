% CORPORATE BOND LIQUIDITY
% Financial Engineering: Politecnico Milano
%
% In order to run the script:
% >> runFinalProject_Group6_1_2
%
% Last Modified: 14.06.2019
%
clc
clear all
close all
addpath(genpath('.'))

%% Loading data and settings

% modified in order to use the same format date of the spreadsheet
formatData = 'dd/mm/yyyy';

[datesSet_dirty, ratesSet, normal_vols] = ...
  readExcelData('curves20150910_project.xlsx',formatData);


% moving ('modified following') the dates that fall in a non-business day
datesSet=clean_date(datesSet_dirty);

%% Discount bootstraps

[dates_EONIA, discounts_EONIA] = bootstrap_OIS(datesSet, ratesSet);

[dates_pseudo,discounts_pseudo] = ...
   bootstrap_pseudo(datesSet,ratesSet, dates_EONIA, discounts_EONIA);

%% Plots of discount curves

plot_discount(dates_EONIA, discounts_EONIA, dates_pseudo, ...
    discounts_pseudo, formatData);

%% Calibration of volatility parameters

[a,sigma,gamma] = calibrate_model(dates_EONIA, discounts_EONIA, ...
    dates_pseudo, discounts_pseudo, normal_vols );

function [dates, rates, normal_vols] = readExcelData( filename, formatData)

% Read data from excel, returning mid-prices and relevant dates.
% All input rates are in % units
%__________________________________________________________________________
% INPUT
% - filename:       excel file name where data are stored;
% - formatData:     data format in Excel.
%--------------------------------------------------------------------------
% OUTPUT
% - dates:      struct with the following fields:
%                       - settlementDate;
%                       - deposDates;
%                       - futuresDates;
%                       - swapDates;
% - rates:      struct with the following fields:
%                       - deposRates;
%                       - futuresRates;
%                       - swapRates.
%--------------------------------------------------------------------------
% Last Modified: 28.05.2019
%__________________________________________________________________________

%% Dates from Excel

%Settlement date
[~, settlement] = xlsread(filename, 1, 'H5');
%Date conversion
dates.settlement = datenum(settlement, formatData);

%Dates relative to OIS
[~, dates_OIS] = xlsread(filename, 1, 'B2:B19');
dates.OIS = datenum(dates_OIS, formatData);

%Dates relative to futures: calc start & end
[~, date_fra_read] = xlsread(filename, 3, 'C2:D10');
numberFRA = size(date_fra_read,1);

dates.fra=ones(numberFRA,2);
dates.fra(:,1) = datenum(date_fra_read(:,1), formatData);
dates.fra(:,2) = datenum(date_fra_read(:,2), formatData);

% Date relative to Euribor6m
[~, date_eur] = xlsread(filename, 2, 'B2');
dates.euribor= datenum(date_eur, formatData);

%Date relative to swaps: expiry dates
[~, date_swaps] = xlsread(filename, 2, 'B3:B14');
dates.swaps = datenum(date_swaps, formatData);

%% Rates from Excel (Bids & Asks)

%OIS (bid and ask)
tassi_depositi = xlsread(filename, 1, 'E2:E19');
rates.OIS = tassi_depositi / 100;

%FRA (bid and ask)
tassi_fra = xlsread(filename, 3, 'G2:G10');
%Rates from futures
rates.fra = tassi_fra/100;

%Euribor6m
euribor6m = xlsread(filename, 2, 'C2');
rates.Euribor6m = euribor6m/100;

%Swaps
tassi_swaps = xlsread(filename, 2, 'C3:C14');
rates.swaps = tassi_swaps / 100;

%% Normal vols

%diagonal normal vols for swaptions
n_vols = xlsread(filename, 4, 'D3:D11');
normal_vols.values = n_vols/ 10000;
normal_vols.expiry=xlsread(filename, 4, 'B3:B11'); % in years
normal_vols.tenor=xlsread(filename, 4, 'C3:C11'); % in years

end
function [datesSet_2, coupon, clean_price, surv_probs] = buildStruct()

% Create the struct containing data of BNPP and Santander on 10 September
% 2015, according to paper [2].
%__________________________________________________________________________
% INPUT
%   no inputs.
%--------------------------------------------------------------------------
% OUTPUT
% - datesSet_2:     struct containing maturity of liquid bonds;
% - coupon:         struct containing coupons of liquid bonds;
% - clean_price:    struct containing end-of-day mid-prices of bonds;
% - surv_probs:     survival probability.
%--------------------------------------------------------------------------
% Last Modified: 05.06.2019
%__________________________________________________________________________
%
% REMARK: we assume to have a face value equal to 1, so everything is
%         rescaled accordingly to this assumption.

%% BNPP data
% maturities
dates_BNPP=['27-Nov-2017';
            '12-Mar-2018';
            '21-Nov-2018';
            '28-Jan-2019';
            '23-Aug-2019';
            '13-Jan-2021';
            '24-Oct-2022';
            '20-May-2024'];
datesSet_2.BNPP=datenum(dates_BNPP);

%coupons
coupon_BNPP=[2.875; 1.5; 1.375; 2; 2.5; 2.25; 2.875; 2.375];
coupon.BNPP=coupon_BNPP/100; %see REMARK

%clean prices
clean_price_BNPP=[105.575; 
                  102.768; 
                  102.555; 
                  104.536; 
                  106.927;
                  106.083; 
                  110.281; 
                  106.007];
clean_price.BNPP=clean_price_BNPP/100; %see REMARK

% survival probabilities
default_probs=[1.27e-4; 
               1.80e-4]; %for 2w and 2m ttl
surv_probs.BNPP=1-default_probs;

%% Santander data
% dates
dates_Santander=['27-Mar-2017';
                 '04-Oct-2017';
                 '15-Jan-2018';
                 '20-Apr-2018';
                 '14-Jan-2019';
                 '13-Jan-2020';
                 '24-Jan-2020';
                 '14-Jan-2022';
                 '10-Mar-2025'];
datesSet_2.Santander=datenum(dates_Santander);

% coupons
coupon_Santander=[4.000; 4.125; 1.750; 0.625; 2.000; 0.875; 4.000;...
                    1.125; 1.125];
coupon.Santander=coupon_Santander/100; %see REMARK

% clean prices
clean_price_Santander=[105.372; 
                       107.358; 
                       102.766; 
                       99.885; 
                       103.984;
                       99.500; 
                       112.836; 
                       98.166; 
                       93.261];
clean_price.Santander=clean_price_Santander/100; %see REMARK

% survival probabilities
default_probs=[5.42e-4; 
               7.73e-4]; %for 2w and 2m ttl
surv_probs.Santander=1-default_probs;

end
function [coupon_payment_dates, idx1, idx1_cs, idx2] = ...
    create_vector_payment_dates(datesSet,settle,ttl)
            
% Build a unique vector containing all coupon payment dates, and return 
% some indexes in order to handle it.
%__________________________________________________________________________
% INPUT
% - datesSet:       vector containing all the maturity dates for the bonds;
% - settle:         settlement date for the basket of bonds;
% - ttl:            time-to-liquidate.
%--------------------------------------------------------------------------
% OUTPUT
% - coupon_payment_dates:  vector that contains all the (unsorted) coupon
%                          payment dates of the bonds (even the ones
%                          between the settlement date and the ttl);
% - idx1:                  vector that contains the number of coupons 
%                          paid for each bond. 
%                          NB: for practical use, the first element is
%                          set to 0, so the first bond will pay idx1(2)
%                          coupons, the second one will pay idx1(3) coupons
%                          and so on.
% - idx1_cs:               cumulative sum of idx_1. This means that
%                          payment dates of the i-th coupon will be
%                          coupon_payment_dates(idx1_cs(i)+1:idx1_cs(i+1))
%                          (indeed if we have N bonds, length(idx1_cs)=N); 
% - idx2:                  vector that contains the number of coupons paid
%                          between the settlement date and the ttl. The
%                          notation used is the same as idx1 and idx1_cs,
%                          with the first element idx2(1) set to 0.
%                          [In our case, since each bond pays annual coupon
%                          and we consider ttl=2weeks or ttl=2months, idx2
%                          is a binary variable that holds 1 if the coupon 
%                          is paid in between and 0 otherwise].
%--------------------------------------------------------------------------
% Functions used: dateMoveVec.
%--------------------------------------------------------------------------
% Last Modified: 08.06.2019
%__________________________________________________________________________

%% Settings

n_bonds=length(datesSet); %number of bonds
num_coupon_payments=floor(yearfrac(settle,datesSet));

coupon_payment_dates_cell=cell(n_bonds,1); %inizialization
idx1=zeros(n_bonds+1,1); %+1 for smarter use in the future
idx2=zeros(n_bonds+1,1); %+1 to use the same notation

%% Creating the vector and indexes

for k=1:n_bonds
	% we compute coupon payment dates
	coupon_payment_dates_cell{k}= dateMoveVec(datesSet(k), ...
            'y', -num_coupon_payments(k):0, 'MF', eurCalendar);
        
    % we store the number of coupons that will be paid for the i-th bond
	idx1(k+1)=length(coupon_payment_dates_cell{k});
    
	% we want to know how many of these coupons fall between settlement and
	% ttl (in our model, at most one, in the worst cases)
	idx2(k+1)=length( ...
        coupon_payment_dates_cell{k}(coupon_payment_dates_cell{k}<=ttl) );
end

% doing the cumulative sum, as explained in the header
idx1_cs=cumsum(idx1);

% converting the dynamic array in a static array
coupon_payment_dates=cell2mat(coupon_payment_dates_cell);

end
function [clean_datesSet] = clean_date(datesSet)

% Take dates and check if they are business dates: if so, they do not
% change, otherwise they become business dates according to 'Modified
% Following' day-convention. The default calendar is eurCalendar.
%__________________________________________________________________________
% INPUT:
% - datesSet:           initial datesSet, in which some dates may fall in
%                       non-business days. It must contain the following
%                       fields:
%                               -settlement
%                               -OIS
%                               -fra
%                               -Euribor
%                               -swaps
%--------------------------------------------------------------------------
% OUTPUT:
% - clean_datesSet:     final datesSet, where all the dates are business
%                       days.
%--------------------------------------------------------------------------
% Last Modified: 28.05.2019
%__________________________________________________________________________

%% Settings

% creating a structure with the same fields and number of elements
clean_datesSet=datesSet;

%% Cleaning fras

clean_datesSet.fra(isbusday(datesSet.fra,eurCalendar)==0) = ...
    busdate(datesSet.fra( isbusday(datesSet.fra,eurCalendar)==0 ), ...
            'modifiedfollow',eurCalendar); 

%% Cleaning OIS

clean_datesSet.OIS(isbusday(datesSet.OIS,eurCalendar)==0) = ...
    busdate(datesSet.OIS( isbusday(datesSet.OIS,eurCalendar)==0 ), ...
            'modifiedfollow',eurCalendar); 

%% Cleaning swaps

clean_datesSet.swaps(isbusday(datesSet.swaps,eurCalendar)==0) = ...
    busdate(datesSet.swaps( isbusday(datesSet.swaps,eurCalendar)==0 ), ...
            'modifiedfollow',eurCalendar); 

end
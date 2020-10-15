function Dates = tenordates(startdate,Nsteps,datepart,num,...
                            businessdayconvention,market)

% Find relevant subsequent dates from a starting date according to a
% regular timestep between dates and adjusted according to a specified
% business day convention and a calendar of holidays.
%__________________________________________________________________________
% INPUT:
% - startdate:              Start date [integer];
% - Nsteps:                 Number of time units to move forward from the
%                           startdate;
% - datepart:               Time unit of the movement: 
%                           - 'd': day;
%                           - 'w': week;
%                           - 'm': month;
%                           - 'y': year.
% - num:                    Number of time units per each timestep;
% - businessdayconvention:  Business Day Convention: 
%                           - 'F': forward;
%                           - 'P': previous;
%                           - 'MF': modified forward;
%                           - 'MP': modified previous;
%                           - 'U': unmodified.
% - market:                 Market: target vector of holidays.
%--------------------------------------------------------------------------
% OUTPUT:
% - Dates:                    Row vector of subsequent dates from the 
%                             startdate arranged according the input num, 
%                             businessdayconvention and vector of holidays.
%--------------------------------------------------------------------------
% Functions used: dateMoveVec.
%__________________________________________________________________________


% Pre-allocate the vector of the relevant dates to be computed
Dates = nan(Nsteps,1) ;

% Find the relevant dates 
for i=1:Nsteps
    Dates(i) = dateMoveVec(startdate,datepart,num*i,...
                             businessdayconvention,market) ;
end

end
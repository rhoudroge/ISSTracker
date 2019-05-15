function dates = getStatesDates( orbitalData )
% GETSTATESDATES return the dates of the state vectors in the orbital data
%   dates = getStatesDates(orbitalData) returns a cell of AbsoluteDate
%   objects that represents the dates (UTC) of all the state vectors
%   contained in orbitalData

import org.orekit.time.*;

% get UTC time scale
utc = TimeScalesFactory.getUTC;

n = length(orbitalData);
dates = cell(n, 1);
for i = 1:n
    % read ith state vector
    c = orbitalData(i).time;
    dates{i} = AbsoluteDate(c.year, c.month, c.day, c.hour, c.minute, ...
        c.second + c.msecond / 1000, utc);
end

end


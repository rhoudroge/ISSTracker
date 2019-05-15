
function [position, velocity, initialDate] = getStateVector(timeNow)
% GETSTATEVECTOR Returns the state vector.
%
% [pos, vel, time] = getStateVector(time) returns the state vector
% (including position, velocity and time) on the same day as input time

import java.lang.Integer;
import java.lang.Math;
import java.lang.System;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.TimeZone;

import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.orekit.time.*;
import org.orekit.data.*;

% Get the UTC time scale
utc = TimeScalesFactory.getUTC;

% The Nasa feed
if simulation_parameters.useSavedFile
    fprintf('     Warning: using local saved file (may be old!)');
    nasaData = getNasaSavedFeed(simulation_parameters.savedFileName);
else
    nasaData = getNasaActualFeed;
end

% get orbital data from feed
orbitalData = getNasaOrbitalData(nasaData);

% Look for todays ephemeris
daysDetected = arrayfun(@(x)orbitalData(x).time.day, 1:length(orbitalData));

% Todays ephemeris
dayToday = timeNow.getComponents(utc).getDate.getDay;
temp = orbitalData(daysDetected == dayToday);

if isempty(temp)
    error('No state vector for given day... Exiting.');
else
    temp = temp(end);
    
    initialDate = AbsoluteDate(temp.time.year, temp.time.month, ...
        temp.time.day, temp.time.hour, temp.time.minute,...
        temp.time.second + temp.time.msecond/1000, utc);
    
    % position and velocity of ISS
    position = Vector3D(temp.data.X, temp.data.Y, temp.data.Z);
    velocity = Vector3D(temp.data.XDot, temp.data.YDot, temp.data.ZDot);
    
    fprintf('     Using state vector dated %s\n', char(initialDate));
end

end
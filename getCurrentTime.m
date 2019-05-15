% getCurrentTime MATLAB code for getCurrentTime.m
%      This function calulates the current time in UTC. It can either be
%      simulated time or real time.
%
%      currentTime = getCurrentTime(TimeScale) return the current time,
%      either simulated or real time.
%
%      currentTime = getCurrentTime(TimeScale, ~) returns the real time.
%
%      getCurrentTime called without any arguments clears the internal
%      persistent variable required by simulated time
%
% See also: simulation_parameters
function currentTime = getCurrentTime(varargin)

persistent firstRun;
persistent origin;

if nargin == 1
    
    utc = varargin{1};
    
    if isempty(firstRun)
        origin = getTimeNow(utc);
        firstRun = false;
    end

    timeNow = getTimeNow(utc);
    
    if simulation_parameters.simulatedTime
        shift = timeNow.durationFrom(origin) * simulation_parameters.speedFactor;
        currentTime = getSimulatedTime(shift, utc);
    else
        currentTime = getTimeNow(utc);
    end
    
elseif nargin == 2
    
    currentTime = getTimeNow(varargin{1});
    
else
    clear origin
    clear firstRun
    getSimulatedTime;
end

end

function timeNow = getTimeNow(utc)


% - Java imports
import java.lang.Math;
import java.lang.System;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.TimeZone;
import org.orekit.time.*;

% - Date of now
date = Date();
calendar = GregorianCalendar.getInstance(TimeZone.getTimeZone('UTC')); % creates a new calendar instance
calendar.setTime(date);   % assigns calendar to given date
year = calendar.get(Calendar.YEAR);
month = calendar.get(Calendar.MONTH)+1;
day = calendar.get(Calendar.DAY_OF_MONTH);
hour = calendar.get(Calendar.HOUR_OF_DAY); % gets hour in 24h format
minute = calendar.get(Calendar.MINUTE);
second = calendar.get(Calendar.SECOND) + calendar.get(Calendar.MILLISECOND)/1000 ;

timeNow =  AbsoluteDate(DateComponents(year, month, day),...
    TimeComponents(hour, minute, second), utc).shiftedBy( ...
    3600 * simulation_parameters.timeShift);

end

function simulatedTime = getSimulatedTime(varargin)

import org.orekit.time.*;

persistent simulationStart;

if nargin == 2
    shift = varargin{1};
    utc = varargin{2};
    
    if isempty(simulationStart)
        start = simulation_parameters.simulationStart;
        startDateComponents = DateComponents(start(1), start(2), start(3));
        startTimeComponents = TimeComponents(start(4), start(5), start(6));
        simulationStart = AbsoluteDate(startDateComponents, startTimeComponents, ...
            utc);
    end

    simulatedTime = simulationStart.shiftedBy(shift);
    
else
    clear simulationStart
end

end
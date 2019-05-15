function computePasses(varargin)
% COMPUTEPASSES Generates a passes prediction file
%   computePasses generates a passes prediction file for the next 48 hours
%   computePasses(hours) generates a passes prediction file for the number
%   of hours specified by the user
%   computePasses(date, hours) generates a passes prediction file for
%   the numberof hours specified by the user starting on the given
%   startDate [yyyy mm dd hh mm sss] in UTC
% Author : Rami Houdroge
% Version : 1.0.0
% Created : 2011
% Revision : $Id: computePasses.m 10 2016-05-08 21:42:58Z rami $

%% Display message
fprintf('Computing passes...\n');

% load java libraries
loadLibraries;
configureLibraries;

%% Load Java classes
% From java
import java.lang.Math;
import java.lang.System;
import java.io.File;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.TimeZone;

% From the Apache Commons Math Project
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.apache.commons.math3.util.*;
import org.apache.commons.math3.ode.nonstiff.*;

% From the ORbit Extrapolation KIT
import org.orekit.bodies.*;
import org.orekit.data.*;
import org.orekit.errors.*;
import org.orekit.frames.*;
import org.orekit.forces.*;
import org.orekit.forces.gravity.*;
import org.orekit.forces.gravity.potential.*;
import org.orekit.forces.radiation.*;
import org.orekit.forces.drag.*;
import org.orekit.orbits.*;
import org.orekit.propagation.*;
import org.orekit.propagation.events.*;
import org.orekit.propagation.numerical.*;
import org.orekit.time.*;
import org.orekit.tle.*;
import org.orekit.utils.*;

% customised orekit elements
import org.orekit.orekit_custom.*;

%% Parse input arguments
utc = TimeScalesFactory.getUTC();
if nargin == 1
    hours = varargin{1};
    timeNow = getCurrentTime(utc, 1);
elseif nargin == 2
    hours = varargin{1};
    dv = varargin{2};
    timeNow = AbsoluteDate(dv(1), dv(2), dv(3), dv(4), dv(5), dv(6), utc);
else
    hours = 48;
    timeNow = getCurrentTime(utc, 1);
end

%% Fetch and parse NASA ephemeris
% get the closest state vector
[pos, vel, time] = getStateVector(timeNow);

% get the associated propagator
propagator = getPropagator(pos, vel, time);
propagator.propagate(timeNow);

% create the user point on the globe
lon = deg2rad(simulation_parameters.userLon);
lat = deg2rad(simulation_parameters.userLat);
alt = simulation_parameters.userAlt;
name = simulation_parameters.userLoc;
coordinates = GeodeticPoint(lat, lon, alt);

% create the rotating non intertial frame
itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, true);

% create a representation of Earth
ae = Constants.GRIM5C1_EARTH_EQUATORIAL_RADIUS;
f = Constants.GRIM5C1_EARTH_FLATTENING;
earth = OneAxisEllipsoid(ae, f, itrf);

% create the Topocentric frame associated to the user viewpoint
topoFrame = TopocentricFrame(earth, coordinates, name);

% create an elevation detector
elevationDetector = CustomElevationDetector(topoFrame);
propagator.addEventDetector(elevationDetector);

% propagate
targetDate = timeNow.shiftedBy(3600 * hours);
propagator.propagate(targetDate);

% get list of events
eventTimes = elevationDetector(1).getEventTimes();
eventDir = elevationDetector(1).getEventDirection();

% start and end times
t_0 = char(timeNow(1).toString);
t_1 = char(targetDate(1).toString);

% same times but adapted for file names
t0 = date2char( timeNow(1) );
t1 = date2char( targetDate(1) );

% open target file
filename = fullfile('out', ['passes_' t0 '_' t1 '.txt']);
fprintf('Writing to file : %s\n', filename);
fileID = fopen(filename,'w');

% print header
fprintf(fileID, 'Pass predictions generated    %s UTC\n', t_0);
fprintf(fileID, '  from                        %s UTC\n', t_0);
fprintf(fileID, '  to                          %s UTC\n', t_1);
fprintf(fileID, 'User location                 %s\n', name);
fprintf(fileID, '  coordinates                 %f°, %f°, %f m\n\n', ...
    rad2deg(lon), rad2deg(lat), alt);

% print events
if length(eventTimes) >= 1
    for i=1:length(eventTimes)
        if eventDir(i) == 1
            str = 'PASS START';
        else
            str = 'PASS END';
        end
        fprintf(fileID, '%s       %s\n',char(eventTimes(i)), str);
    end
else
    fprintf(fileID, 'NO PASSES OVER SELECTED PERIOD');
end

% close file
fclose(fileID);

fprintf('Done...\n');

end

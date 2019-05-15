
function computeDeviations
%
% Provides 24h propagation errors for each bulletin in the NASA feed.
% Author : Rami Houdroge
% Version : 1.0.0
% Created : 2011
% Revision : $Id: computeDeviations.m 10 2016-05-08 21:42:58Z rami $

clc
loadLibraries;

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
import org.orekit.propagation.numerical.*;
import org.orekit.time.*;
import org.orekit.tle.*;
import org.orekit.utils.*;

configureLibraries;

data = getNasaOrbitalData(getNasaActualFeed);

utc = TimeScalesFactory.getUTC();

n = getCurrentTime(utc, 1);

n1 = date2char(n);
filename = fullfile('out', ...
    ['errors_', n1, '.txt']);

fid = fopen(filename, 'w');
fprintf('Writing to file : %s\n\n', filename);
fprintf(fid, 'Deviations generated %s UTC\n\n', char(n.toString));

for start = 1:length(data)-1
    
    clear propagator
    
    pt1 = data(start);
    pt2 = data(start + 1);
    
    pt1pv = pt1.data;
    pt2pv = pt2.data;
    
    pt1t = pt1.time;
    pt2t = pt2.time;
    
    p1 = Vector3D(pt1pv.X, pt1pv.Y, pt1pv.Z);
    v1 = Vector3D(pt1pv.XDot, pt1pv.YDot, pt1pv.ZDot);
    
    p2 = Vector3D(pt2pv.X, pt2pv.Y, pt2pv.Z);
    v2 = Vector3D(pt2pv.XDot, pt2pv.YDot, pt2pv.ZDot);
    
    d1 = AbsoluteDate(DateComponents(pt1t.year, pt1t.month, pt1t.day), TimeComponents(pt1t.hour, pt1t.minute, pt1t.second), utc);
    d2 = AbsoluteDate(DateComponents(pt2t.year, pt2t.month, pt2t.day), TimeComponents(pt2t.hour, pt2t.minute, pt2t.second), utc);
    
    propagator = getPropagator(p1, v1, d1);

    %% propagate
    startT = System.currentTimeMillis;
    finalSpacecraftState = propagator.propagate(d2);
    endT = System.currentTimeMillis;
    
    pv = finalSpacecraftState.getPVCoordinates();
    p = pv.getPosition();
    v = pv.getVelocity();
    
    ep  = p.subtract(p2).getNorm();
    ev = v.subtract(v2).getNorm();
    
    fprintf(fid, 'Run #%i : %.2f s\n', start, (endT - startT)/1000);
    fprintf(fid, 'Orbit bulletin dated %s UTC\n', char(d1.toString()));
    fprintf(fid, 'Errors after 24h : %d m and %d m/s\n\n', ep, ev);
end

fclose(fid);

end
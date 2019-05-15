% ISSTracker simulation parameters. User can customize the simulation to
% display longer or shorter ground track, use simulated time or real time,
% change his/her location on the surface and display the coast lines or
% Earth picture from the Blue Marble collection.
%
% Author : Rami Houdroge
% Version : 1.1
% Created : 2011
% Revision : $Id: simulation_parameters.m 10 2016-05-08 21:42:58Z rami $
classdef simulation_parameters
    properties (Constant)
        % Nasa ephemeris
        useSavedFile = false;
        savedFileName = 'nasa_feed_20151014.txt';
        url = 'https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html';
                
        % Use real time or simulated time (UTC)
        simulatedTime = true;
        simulationStart = [2019 5 15 21 48 10];
        speedFactor = 1;
        
        % In case we use real time, we can shift the date (in hours)
        timeShift = 0;
        
        % ephemeris times - ephemeris will be valid between
        % now - eph start hours
        % now + eph end hours
        ephStart = -2;
        ephEnd = 5;
        
        % ground trace points - ground trace to be shown between
        % now - gtStart hours
        % now + gtEnd hours
        %
        % make sure that    ephStart < gtStart
        %       and that      ephEnd > gtEnd
        points = 750;
        gtStart = -.5;
        gtEnd = 4.5;
        
        % footprint points
        footPrint = 50;
        
        % user location - degrees
        % Toulouse
        %         userLat = -20;
        %         userLon = 130;
        %         userAlt = 52;
        
        % Darmstadt
        userLoc = 'Darmstadt, DE';
        userLat = 46.8667;
        userLon =  8.65;
        userAlt = 0;

        % propagator
        minStep  = 0.1;
        maxStep  = 50;
        absTolerance = [1e-10 1e-10 1e-10 1e-13 1e-13 1e-13 1e-10];
        relTolerance = [1e-8 1e-8 1e-8 1e-10 1e-10 1e-10 1e-7];
        
        % 3D viewpoint
        latCoef = .5;
        
        % Render Earth
        renderEarth = true;
        
        % Default 3D plot : earth or sky
        default = 'sky';

        % sky map catalog file name and min magnitude
        skyMapCatalog = 'sao6.csv';
        catalogRegEx = '%d;%f;%f;%f;%f;%f';
        minMag = 4;
        
        % number of tick marks to display form each pass on the sky map
        tickMarks = 2;
    end
end


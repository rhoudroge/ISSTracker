%%
% ISSTracker physical parameters
% Author : Rami Houdroge
% Version : 1.1
% Created : 2011
% Revision : $Id: physical_parameters.m 10 2016-05-08 21:42:58Z rami $
%%
classdef physical_parameters
    properties (Constant)
        % earth
        f = 1/298.257223563;
        ae = 6378136.46;
        mu = 3.986004415e14;
        
        % degree and order of earth potential
        potentialFlag = true;
        degree = 10;
        order = 10;
        
        % third bodies
        sunAttraction = true;
        moonAttraction = true;
        innerBodies = false;
        outerBodies = false;
        
        % drag
        dragFlag = false;
        solarActivityFileName = 'Oct2014F10.txt';
        
        
        % solar radiation pressure
        srpFlag = false;

        % ISS
        pressurizedVolume = 937; % m^3
        dim = physical_parameters.pressurizedVolume^(1/3);
        solarArrayArea = 34 * 12 * 8;
        solarArrayRotAxis = [0, 0, 1];
        dragCoefficient = 2;
        absorptionCoefficient = .1;
        reflectionCoefficient = .9;
        

    end
    
end

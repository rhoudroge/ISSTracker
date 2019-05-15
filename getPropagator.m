%%
% ISSTracker propagator definition. File separate from ISSTracker.m for
% test purposes.
% Author : Rami Houdroge
% Version : 1.1
% Created : 2011
% Revision : $Id: getPropagator.m 10 2016-05-08 21:42:58Z rami $
%%
% --- Initialization of Orekit Data / Executes on load
function propagator = getPropagator(p, v, iD)
% s = 'call'
% - Java imports

% From java
import java.lang.Math;
import java.lang.System;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.TimeZone;
import java.io.File;

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
import org.orekit.forces.drag.MarshallSolarActivityFutureEstimation;
import org.orekit.orbits.*;
import org.orekit.propagation.*;
import org.orekit.propagation.numerical.*;
import org.orekit.time.*;
import org.orekit.tle.*;
import org.orekit.utils.*;

% data pointer
% System.setProperty(DataProvidersManager.OREKIT_DATA_PATH, 'data');
DM=DataProvidersManager.getInstance();
crawler=DirectoryCrawler(File(fullfile(cd, 'data')));
DM.clearProviders();
DM.addProvider(crawler);

% initialize UTC time scale
utc = TimeScalesFactory.getUTC();
orekitData.utc = utc;

% - GUI Data pointer
hObject = findall(0,'Tag','ISSTracker');
if ~isempty(hObject)
    handles = guidata(hObject);
end

% - ISS live data feed

% - Frames

% Frames used (EME2000 and ITRF)
orekitData.frames.EME2000 = FramesFactory.getEME2000();
orekitData.frames.ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, false);

% Earth Station Location (Toulouse)
lat = simulation_parameters.userLat;
lon = simulation_parameters.userLon;
alt = simulation_parameters.userAlt;

stationPoint = GeodeticPoint(FastMath.toRadians(lat),...
    FastMath.toRadians(lon), alt);


% terrestrial frame
orekitData.earthBody = OneAxisEllipsoid(physical_parameters.ae,...
    physical_parameters.f, orekitData.frames.ITRF);
orekitData.frames.stationFrame = TopocentricFrame(orekitData.earthBody,...
    stationPoint, 'Toulouse');

% Initial Orbit
initialOrbit = CartesianOrbit(PVCoordinates(p, v),...
    orekitData(1).frames.EME2000, iD, physical_parameters.mu);
initialState = SpacecraftState(initialOrbit);


% Variable Step Dormand Prince Integration
integrator = DormandPrince853Integrator(simulation_parameters.minStep,...
    simulation_parameters.maxStep, simulation_parameters.absTolerance, simulation_parameters.relTolerance);
propagator = NumericalPropagator(integrator);
propagator.setInitialState(initialState);
propagator.setEphemerisMode();


% - Forces applied on ISS
% Earth Gravity - Holmes Featherstone Attraction Model with EIGEN-5C Gravity Model

% Gravitational Forces
if physical_parameters.potentialFlag
    % get the coefficients of the chosen model
    GravityFieldFactory.clearPotentialCoefficientsReaders;
    coefs = GravityFieldFactory.getNormalizedProvider(physical_parameters.degree, physical_parameters.order);
    % create the force model instance
    earG = HolmesFeatherstoneAttractionModel(orekitData(1).frames.ITRF, coefs);
    % add to propagator
    propagator.addForceModel(earG);
end

% sun model
sun = CelestialBodyFactory.getSun();

% ISS geometric model
solarArrayRotAxis = Vector3D(physical_parameters.solarArrayRotAxis);
dim = physical_parameters.dim;
ISSmodel = BoxAndSolarArraySpacecraft(dim, dim, dim, ...
    sun, physical_parameters.solarArrayArea, solarArrayRotAxis, ...
    physical_parameters.dragCoefficient, physical_parameters.absorptionCoefficient, physical_parameters.reflectionCoefficient);

% Atmosphere drag force and sun radiation pressure
% Both forces require a spacecraft model. A box and solar array
% approximation is considered here with the same physical parameters as the
if physical_parameters.dragFlag
      
    % solar activity strength level
    enumLoc = 'org.orekit.forces.drag.MarshallSolarActivityFutureEstimation$StrengthLevel';
    enumValues = cell(javaMethod('values', enumLoc));
    average = enumValues{2};

    % atmosphere model
    msafe = MarshallSolarActivityFutureEstimation(physical_parameters.solarActivityFileName, ...
        average);
    DataProvidersManager.getInstance().feed(msafe.getSupportedNames(), msafe);
    atmosphere = DTM2000(msafe, sun, orekitData(1).earthBody);
    
    dragForce = DragForce(atmosphere, ISSmodel);
    propagator.addForceModel(dragForce);
end

% Solar radiation pressure
if physical_parameters.srpFlag
    sunRadiation = SolarRadiationPressure(sun, 1.392e6, ISSmodel);
    propagator.addForceModel(sunRadiation);
end



% - Add forces to propagator

if physical_parameters.sunAttraction
    sunG = ThirdBodyAttraction(CelestialBodyFactory.getSun());
    propagator.addForceModel(sunG);
end
if physical_parameters.moonAttraction
    mooG = ThirdBodyAttraction(CelestialBodyFactory.getMoon());
    propagator.addForceModel(mooG);
end
if physical_parameters.innerBodies
    merG = ThirdBodyAttraction(CelestialBodyFactory.getMercury());
    venG = ThirdBodyAttraction(CelestialBodyFactory.getVenus());
    marG = ThirdBodyAttraction(CelestialBodyFactory.getMars());
    propagator.addForceModel(merG);
    propagator.addForceModel(venG);
    propagator.addForceModel(marG);
end
if physical_parameters.outerBodies
    jupG = ThirdBodyAttraction(CelestialBodyFactory.getJupiter());
    satG = ThirdBodyAttraction(CelestialBodyFactory.getSaturn());
    uraG = ThirdBodyAttraction(CelestialBodyFactory.getUranus());
    nepG = ThirdBodyAttraction(CelestialBodyFactory.getNeptune());
    pluG = ThirdBodyAttraction(CelestialBodyFactory.getPluto());
    propagator.addForceModel(jupG);
    propagator.addForceModel(satG);
    propagator.addForceModel(uraG);
    propagator.addForceModel(nepG);
    propagator.addForceModel(pluG);
end

propagator.setOrbitType(OrbitType.CARTESIAN);

% - Save data structure in GUI
set(hObject, 'UserData', orekitData);

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                                   :          ....
%                                .. .~  ..:77..$II$I+.
%                                 :Z=   .$I7ZII+.,II=7?7,.
%                            .     ,=.   ..,?I?7II:..777$7+$.
%                           ,N      .~     ...$II+~I?,.~=7$O+.
%       .  ~77~I77.         .=,.,.   +~=..    ...7II,Z....7+I77II7.
% =77.$77. .O7II7I7..      ..$$$N, .,=N:7~:.     ...+,I$7I7..Z$7IZ$7O.
% .777.777.  ?777$777.       77Z=~..$+7~$I~N....77~D8.?7?+7I7,..Z777Z$IZ..
% .$I$IZ77I   ~77I.77I..    .7$$?7M7ZZ$~Z.+~?: .... 8,..77II$7I7...Z77$$$7O..
%  .7$7.7777. .~777.7I77.. ...?Z~I8=~~~~Z 7O    =:$..O?. .$I77?7$7=..,I$7$Z$7O..
%  ..$$7 777$.  .I77.7$~...I.~:?I$~.I$O:.     ..,888=8M    ..777$7I77...II$$Z$$$7.
%    $$II.I777  ..:$=+?+~7:~8O ,==..~:M+,.+N+~+?77~8,.       ..$$$7I77$Z...7$,..
%    .777$.7?...~,=77I7O7I77..+7N.:II,,,:77D8~..                .~777777ID.
%     ..O$7II+$?=ZN.7I$7.777I...$$$$877I::DN8                    ...$D..
%     .+7I7=$777?. ..7$$I.77I7..Z...Z+7D..
%     ..$$$7.Z777+.  .7777.I7?++ ...7D77..  
%       .$777.$777,   .77I$.ZI?+:....M,
%       .?$777,$77I.. .,?III,.?$~+.
%        .$777$.I7$7,.  .III+?.7I=+.
%        ..77$7=.77$7.   .7?Z7=.+II7~.
%          .$7$7..77$$,  ..$$+I=.$+$?I.
%          .$7$77.~7$7$.   .7$?$?.O$O.
%            $777$.O777$..   $7$8.
%            .7$777.7II$7.
%             :$I7I~.$Z...
%             .$7$O..
%
%
% ISSTracker main script
% Author : Rami Houdroge
% Version : 1.1
% Created : 2011
% Revision : $Id: ISSTracker.m 10 2016-05-08 21:42:58Z rami $
%
%% GUIDE methods

% ---
function varargout = ISSTracker(varargin)
% ISSTRACKER MATLAB code for ISSTracker.fig
%      ISSTRACKER, by itself, creates a new ISSTRACKER or raises the existing
%      singleton*.
%
%      H = ISSTRACKER returns the handle to a new ISSTRACKER or the handle to
%      the existing singleton*.
%
%      ISSTRACKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ISSTRACKER.M with the given input arguments.
%
%      ISSTRACKER('Property','Value',...) creates a new ISSTRACKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ISSTracker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ISSTracker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ISSTracker

% Last Modified by GUIDE v2.5 09-Feb-2015 21:58:16

% check that commons math3 and orekit are loaded

if checkLibrariesStatus == 0
    % if not interrupt
    loadLibraries;
    h = warndlg('Loaded required libraries.');
    set(h, 'Tag', 'ISSTW');
    button = findall(0, 'Tag', 'OKButton');
    set(button, 'String', 'Run ISSTracker...');
    pos = get(button, 'Position');
    set(button, 'Position', [pos(1)/1.35 pos(2) pos(3)*2.7 pos(4)]);
    set(button, 'Callback', @warndlgCallback);
    
else
    % otherwise continue
    
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
        'gui_Singleton',  gui_Singleton, ...
        'gui_OpeningFcn', @ISSTracker_OpeningFcn, ...
        'gui_OutputFcn',  @ISSTracker_OutputFcn, ...
        'gui_LayoutFcn',  [] , ...
        'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end
    
    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
end
end

% --- callback for warning dialog
function warndlgCallback(varargin)
close(findall(0, 'Tag', 'ISSTW'));
ISSTracker;
end

% ---
function ISSTracker_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ISSTracker (see VARARGIN)

% Choose default command line output for ISSTracker
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Call main
main(hObject, handles)

% UIWAIT makes ISSTracker wait for user response (see UIRESUME)
% uiwait(handles.ISSTracker);
end

% ---
function varargout = ISSTracker_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

%% Main

% --- Calls data & GUI initialization methods

function main(hObject, handles)

% clear console
clc

% import some Orekit classes
import org.orekit.time.*;
configureLibraries;

% clear persistent variables in getCurrentTime (done by calling function
% w/o arguments)
getCurrentTime
utc = TimeScalesFactory.getUTC;

fprintf('  Launching ISSTracker\n  ====================\n\n')

if simulation_parameters.simulatedTime
    s = simulation_parameters.simulationStart;
    startDate = AbsoluteDate(DateComponents(s(1), s(2), s(3)),...
        TimeComponents(s(4), s(5), s(6)), utc).toString;
    fprintf('   Using simulated time starting : %s\n', char(startDate));
    fprintf('                     speed fator : %f\n\n',...
        char(simulation_parameters.speedFactor));
else
    fprintf('   Using real time \n\n');
end

timeNow = getCurrentTime(utc);

status = get(hObject, 'UserData');
if isempty(status)
    
    set(handles.info, 'UserData', 0);
    set(handles.startTrackerButton, 'UserData', 0);
    
    %% Retrieve and parse online data
    fprintf('   Retrieving ISS Orbital Data...\n');
    [p, v, iD] = getStateVector(timeNow);
    
    %% Propagate orbit bulletin
    propagator = getPropagator(p, v, iD);
    
    % and store generated ephemeris
    orekitData = get(hObject, 'UserData');
    ephemeris = propagate(propagator, orekitData(1).utc);
    orekitData.generatedEphemeris = ephemeris;
    
    
    %% 3D view initialization
    
    % what 3D plot ? 'sky' or 'earth'
    set(handles.uipanel1, 'UserData', simulation_parameters.default);
    
    % Date of now
    currentTime = getCurrentTime(orekitData(1).utc);
    orekitData.startTime = currentTime;
    set(hObject, 'UserData', orekitData);
    
    % plot data
    [~, lon, lat, azs, els, ~, x, y, z, ~, ~, ~, ~, ~, ~, id, times] = ...
        getTimeSeries (currentTime, simulation_parameters.gtStart, simulation_parameters.gtEnd, orekitData, orekitData(1).frames.EME2000);
    % earth rotation
    zcomp = orekitData.frames.ITRF.getTransformTo(orekitData.frames.EME2000, currentTime).getRotation.getAxis.getZ;
    era = orekitData.frames.ITRF.getTransformTo(orekitData.frames.EME2000, currentTime).getRotation.getAngle;
    
    % init method
    if strcmpi(simulation_parameters.default, 'earth')
        earth3D = 'on';
        sky3D = 'off';
    else
        earth3D = 'off';
        sky3D = 'on';
    end
    
    fprintf('   Initializing Sky Map...\n')
    initializeSkyMap(sky3D);
    initializeSkyMapTrail(azs, els, times);
    updateSkyMapISS(azs(id), els(id));
    
    fprintf('   Initializing Earth View...\n')
    initializeEarthView(x, y, z, id, lon, lat, era, zcomp, earth3D);
    
    %% 2D view initialization
    fprintf('   Initializing Ground Track View...\n')
    % plot data
    [t, lon, lat, az, el, alt, ~, ~, ~, a, e, i, pa, ma, raan, id, ~] = ...
        getTimeSeries (currentTime, simulation_parameters.gtStart, simulation_parameters.gtEnd, orekitData, orekitData(1).frames.ITRF);
    
    % init method
    r = sqrt(x(id)^2 + y(id)^2 + z(id)^2);
    initializeGroundTrackView(lon, lat, id, r);
    
    %% Station viewpoint initialization
    % uses same time span as earlier
    fprintf('   Initializing Other Plots...\n\n')
    initializeOtherPlots(az, el, t, a, e, i, pa, ma, raan, alt, id);
    
    fprintf('   Done! Press "Start" to start tracking...\n')
end



end

function status = checkLibrariesStatus()

loaded = javaclasspath('-dynamic');
cmFlag = 0;
orekitFlag = 0;

for k=1:length(loaded)
    if ~isempty(strfind(loaded{k}, 'math'))
        cmFlag = 1;
    end
    if ~isempty(strfind(loaded{k}, 'orekit'))
        orekitFlag = 1;
    end
end
status = cmFlag && orekitFlag;

end

% --- Executes on button press in info.
function info_Callback(~, ~, handles)
% hObject    handle to info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.info, 'UserData', 1);
handleInfoEvent();

end

% --- Handles launching the info gui
function handleInfoEvent()


hObject = findall(0, 'Tag', 'ISSTracker');
handles = guidata(hObject);
wasOn = get(handles.startTrackerButton, 'UserData');
if ~wasOn
    % Otherwise the thread calls this method
    displayInfoScreen();
end

end

% --- Loads information screen
function displayInfoScreen()

hObject = findall(0, 'Tag', 'ISSTracker');
handles = guidata(hObject);
infoFig = findall(0, 'Name', 'ISSTrackerInfo');
orekitData = get(hObject, 'UserData');

% uiwait(handles.ISSTracker)

if isempty(infoFig)
    % size
    bgColor = get(hObject, 'Color');
    fgColor = get(handles.uipanel1, 'ForegroundColor');
    s = get(0, 'ScreenSize');
    w = 670;
    h = 300;
    ox = (s(3)-w)/2;
    oy = 2*(s(4)-h)/3;
    
    % new figure
    infoFig = figure('Position',[ox oy w h], 'Color',...
        bgColor, 'Name', 'ISSTrackerInfo',...
        'NumberTitle', 'off', 'Toolbar', 'none', 'MenuBar', 'none', ...
        'Units', 'normalized', 'Visible', 'off');
    
    % controls
    uicontrol('Style', 'Text', 'String', 'Information & Acknowledgments',...
        'units', 'normalized', 'position', [.2 .80 .6 .1], 'FontWeight', ...
        'bold', 'fontsize', 12, 'Backgroundcolor', bgColor, 'ForegroundColor', ...
        fgColor);
    
    infoStr = ['The ISS Real Time Tracker is a Matlab interface developed ', ...
        'with GUIDE. It allows real time station tracking by retrieving the ', ...
        'ISS orbital ephemeris from the NASA online JPL database and propagating', ...
        ' it using the ORbit Extrapolation KIT', ...
        ' (OREKIT) which is a free and open source low-level space dynamics', ...
        ' library. Many thanks go out to all the teams that have dedicated', ...
        ' much time and effort in order to provide these public data sets and', ...
        ' libraries, making smaller projects like these possible.'];
    
    uicontrol('Style', 'Text', 'String', infoStr,...
        'units', 'normalized', 'position', [.05 .40 .9 .2], ...
        'Backgroundcolor', bgColor, 'ForegroundColor', ...
        fgColor, 'HorizontalAlignment', 'left');
    
    figure(infoFig);
else
    figure(infoFig);
end
set(handles.info, 'UserData', 0);

end

%% GUI initialization and ISS data recovery methods

% --- Propagate orbit bulletin
function ephemeris = propagate(propagator, utc)

fprintf('   Calculating ephemeris...\n')
currentTime = getCurrentTime(utc);

propagator.propagate(currentTime.shiftedBy(simulation_parameters.ephStart * 3600));
propagator.propagate(currentTime.shiftedBy(simulation_parameters.ephEnd * 3600));

ephemeris = propagator.getGeneratedEphemeris();

end

% --- Computes ground trace
function [t, lons, lats, az, el, alt, x, y, z, a, e, i, pa, ma, raan, id, groundDate]...
    = getTimeSeries (time, start, endH, orekitData, frame)

import java.lang.Math;
import java.lang.System;
import org.orekit.time.*;

% time axis
p = simulation_parameters.points;
t = unique(sort([linspace(60 * 60 * start, 60 * 60 * endH, p), 0]));
groundDate = arrayfun(@(x)time.shiftedBy(x), t, 'Un', 0);

% get orbital info for each step
ISSData = cellfun(@(x)getISSData(x, orekitData, frame), groundDate);

% latitude, longitude, azimuth and elevation
lons = [ISSData.lon];
lats = [ISSData.lat];
az = [ISSData.azimuth];
el = [ISSData.elevation];
alt = [ISSData.alt];

% cartesian elements
x = [ISSData.pX];
y = [ISSData.pY];
z = [ISSData.pZ];

% keplerian elements
a = [ISSData.a];
e = [ISSData.e];
i = [ISSData.i];
pa = [ISSData.pa];
ma = [ISSData.ma];
raan = [ISSData.raan];

% start index
id = t == 0;

end

%% Simulation methods

% --- Main thread
function startTracking(hObject, handles)

orekitData = get(hObject, 'UserData');

runStatus = 1;
rtn = 1;
infoRq = 0;

while runStatus && rtn && ~infoRq
    
    % Main "thread" that is running as long as the startTrackerButton
    % is on ('UserData' is 1)
    handles = guidata(hObject);
    runStatus = get(handles.startTrackerButton, 'UserData');
    infoRq = get(handles.info, 'UserData');
    
    % The following two lines allow the user to interact with the GUI
    % during runtime
    pause(.01)
    drawnow;
    
    if runStatus
        
        % - Get updated ISS data
        currentTime = getCurrentTime(orekitData.utc);
        rtn = refreshISS(orekitData, currentTime);
        
    end
    
    % If information screen requested by user, stop running and call
    % appropriate method. Thread is interrupted by ~infoRq = 0
    if infoRq
        set(handles.startTrackerButton, 'UserData', 0);
        set(handles.startTrackerButton, 'String', 'Start');
        set(handles.info, 'UserData', 0);
        
        displayInfoScreen();
    end
    
end


end

% --- Computation of new ISS data
function rtn = refreshISS(orekitData, time, whichFrame)
persistent oldPxs;
if isempty(oldPxs)
    oldPxs = 0;
end

persistent issSky;
if isempty(issSky)
    issSky = 0;
end

persistent sky;
if isempty(sky)
    sky = 0;
end

hObject = findall(0,'Tag','ISSTracker');
if ishghandle(hObject)
    
    % get GUI data
    handles = guidata(hObject);
    
    plot3D = get(handles.uipanel1, 'UserData');
    
    switch plot3D
        case 'earth'
            % earth rotation angle
            zcomp = orekitData.frames.ITRF.getTransformTo(orekitData.frames.EME2000, time).getRotation.getAxis.getZ;
            ERA = orekitData.frames.ITRF.getTransformTo(orekitData.frames.EME2000, time).getRotation.getAngle;
            
            % - Get orbital data for now
            ISSData = getISSData(time, orekitData, orekitData(1).frames.EME2000);
            
            % - Update position in 3D graph
            pos = [ISSData.pX; ISSData.pY; ISSData.pZ];
            s = .05;
            plotX = [0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*s -.025 + pos(1);
            plotY = [0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*s -.025 + pos(2);
            plotZ = [0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*s -.025 + pos(3);
            
            if get(handles.boundRotation, 'Value')
                set(handles.ThDV, 'View', [ISSData.lon+90 + sign(zcomp) * ERA * 180 / pi, ISSData.lat * simulation_parameters.latCoef]);
            end
            
            for i=1:6
                set(handles.cube(i),'XData',plotX(:,i));
                set(handles.cube(i),'YData',plotY(:,i));
                set(handles.cube(i),'ZData',plotZ(:,i));
            end
            
            if simulation_parameters.renderEarth
                
                % Get image width in pixels
                width = size(handles.im, 2);
                
                % determine the angle per pixel
                anglePerPx = 2 * pi / width;
                
                % determine how many pixels need to be moved
                pxs = floor(ERA / anglePerPx);
                
                % - Earth rotation!!
                % Instead of drawing a new sphere each time, we rotate the image
                % (CData)! It is much less intensive computationally
                if (pxs ~= oldPxs);
                    
                    for k=1:3
                        current = handles.im(:,:,k);
                        if zcomp > 0
                            newCdata(:,:,k) = [current(:, width - pxs:width), current(:, 1:width - pxs - 1)];
                        else
                            newCdata(:,:,k) = [current(:, pxs + 1:width), current(:, 1:pxs)];
                        end
                    end
                    set(handles.myEarth, 'CData', newCdata);
                    
                    oldPxs = pxs;
                end
                
                
            else
                theta = handles.coastLineTheta + zcomp * ERA;
                phi = handles.coastLinePhi;
                coastLine = [cos(theta).*cos(phi),...
                    sin(theta).*cos(phi),...
                    -sin(phi)];
                
                set(handles.myEarth,'XData',coastLine(:,1));
                set(handles.myEarth,'YData',coastLine(:,2));
                set(handles.myEarth,'ZData',-coastLine(:,3));
                
            end
            
        case 'sky'
            ISSData = getISSData(time, orekitData, orekitData(1).frames.ITRF);
            
            % - update sky map
            issSky = issSky + 1;
            if issSky == 10
                updateSkyMapISS(ISSData.azimuth, ISSData.elevation);
                issSky = 0;
            end
            
            % - update sky map
            sky = sky + 1;
            if sky == 500
                updateSkyMapStars(orekitData(1).utc);
                sky = 0;
            end
    end
    
    % - Update ground trace
    set(handles.myISS2D, 'XData', ISSData.lon+180);
    set(handles.myISS2D, 'YData', -ISSData.lat);
    
    % - Update footprint
    beta = acos(1 / ISSData.r);
    footLons = ISSData.lon+180 + rad2deg(beta * cos(handles.footprintArray));
    footLats = -ISSData.lat + rad2deg(beta * sin(handles.footprintArray));
    idxs = footLons > 360 | footLons < 0;
    
    ffootLons = [];
    ffootLats = [];
    
    if sum(idxs) > 0
        
        if sum(footLons > 360) >= 1
            % right side
            ffootLons = footLons(idxs) - 360;
            n = length(ffootLons);
            ffootLons = [ffootLons(1:n/2) 0 0 ffootLons(n/2+1:n)];
            ll = rad2deg(beta) * sin(acos((180 - ISSData.lon)/rad2deg(beta)));
            ll1 = -ll - ISSData.lat;
            ll2 =  ll - ISSData.lat;
            footLons(idxs) = 360;
            ffootLats = footLats(idxs);
            ffootLats = [ffootLats(1:n/2) ll2 ll1 ffootLats(n/2+1:n)];
        else
            % left side
            ffootLons = [360 footLons(idxs) + 360 360];
            ll = rad2deg(beta) * sin(acos((ISSData.lon + 180)/rad2deg(beta)));
            ll1 = -ll - ISSData.lat;
            ll2 =  ll - ISSData.lat;
            footLons(idxs) = 0;
            ffootLats = [ll2 footLats(idxs) ll1];
        end
        
    end
    
    set(handles.ISSFootPrint, 'XData', footLons);
    set(handles.ISSFootPrint, 'YData', footLats);
    if isempty(ffootLons)
        set(handles.ISSFootPrintMinor, 'Visible', 'off')
    else
        set(handles.ISSFootPrintMinor, 'Visible', 'on')
    end
    set(handles.ISSFootPrintMinor, 'XData', ffootLons);
    set(handles.ISSFootPrintMinor, 'YData', ffootLats);
    
    % - Update Text Fields
    set(handles.timeInfo, 'String', ISSData.datestr);
    set(handles.a, 'String', sprintf('%.4f', ISSData.a));
    set(handles.e, 'String', sprintf('%.7f', ISSData.e));
    set(handles.i, 'String', sprintf('%.4f', ISSData.i));
    set(handles.pa, 'String', sprintf('%.4f', ISSData.pa));
    set(handles.raan, 'String', sprintf('%.4f', ISSData.raan));
    set(handles.ma, 'String', sprintf('%.4f', mod(ISSData.ma,360)));
    set(handles.lat, 'String', sprintf('%.4f', ISSData.lat));
    set(handles.lon, 'String', sprintf('%.4f', ISSData.lon));
    set(handles.alt, 'String', sprintf('%.4f', ISSData.alt));
    set(handles.az, 'String', sprintf('%.4f', ISSData.azimuth));
    set(handles.el, 'String', sprintf('%.4f', ISSData.elevation));
    set(handles.ra, 'String', sprintf('%.4f', ISSData.range));
    set(handles.vel, 'String', sprintf('%.6f', ISSData.v.getNorm()));
    
    % Switch curves at user request
    axes(handles.TP)
    type = get(handles.whatplot, 'UserData');
    
    current_t = time.durationFrom(orekitData.startTime);
    switch type
        case 'altt'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.alt);
        case 'at'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.a);
        case 'et'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.e);
        case 'az'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.az);
        case 'el'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.el);
        case 'it'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.i);
        case 'pat'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.pa);
        case 'mat'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.ma);
        case 'rat'
            set(handles.ISSSta, 'XData', current_t);
            set(handles.ISSSta, 'YData', ISSData.raan);
        case 'azel'
            set(handles.ISSSta, 'XData', ISSData.azimuth);
            set(handles.ISSSta, 'YData', ISSData.elevation);
    end
    
    
    %     sprintf('az : %.2f, el : %.2f, ra : %.2f', ISSData.azimuth, ISSData.elevation, ISSData.range)
    rtn = 1;
else
    rtn = 0;
    infoFig = findall(0, 'Name', 'ISSTrackerInfo');
    if ~isempty(infoFig)
        figure(infoFig)
    end
end

end

% --- Computation of new ISS data for given date
function ISSData = getISSData(currentDate, orekitData, frame)

% Imports
import java.lang.Math;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.apache.commons.math3.util.*;
import org.orekit.bodies.*;
import org.orekit.data.*;
import org.orekit.errors.*;
import org.orekit.frames.*;
import org.orekit.orbits.*;
import org.orekit.propagation.analytical.*;
import org.orekit.time.*;
import org.orekit.tle.*;
import org.orekit.utils.*;

% Position and Velocity in specified frame
result = orekitData.generatedEphemeris.propagate(currentDate);
finalOrbit = KeplerianOrbit(result.getOrbit());
ISSPV = result.getPVCoordinates(frame);
ISSPosition = ISSPV.getPosition();
ISSVelocity = ISSPV.getVelocity();
ISSLatLon = orekitData.earthBody.transform(ISSPosition,...
    frame, currentDate);


% Get and Store data in ISSData structure
ISSData.datestr = char(currentDate.toString());
ISSData.lat = Math.toDegrees(ISSLatLon.getLatitude());
ISSData.lon = Math.toDegrees(ISSLatLon.getLongitude());
ISSData.alt = ISSLatLon.getAltitude() / 1000;
ISSData.p = ISSPosition;
ISSData.pX = ISSPosition.getX() / physical_parameters.ae;
ISSData.pY = ISSPosition.getY() / physical_parameters.ae;
ISSData.pZ = ISSPosition.getZ() / physical_parameters.ae;
ISSData.r = ISSPosition.getNorm() / physical_parameters.ae;
ISSData.v = ISSVelocity;
ISSData.vX = ISSVelocity.getX();
ISSData.vY = ISSVelocity.getY();
ISSData.vZ = ISSVelocity.getZ();
ISSData.a = finalOrbit.getA();
ISSData.e = finalOrbit.getE();
ISSData.i = Math.toDegrees(finalOrbit.getI());
ISSData.pa = Math.toDegrees(finalOrbit.getPerigeeArgument());
ISSData.ma = Math.toDegrees(finalOrbit.getMeanAnomaly());
ISSData.raan = Math.toDegrees(finalOrbit.getRightAscensionOfAscendingNode());
frames = orekitData.frames(1);
station = frames(1).stationFrame;

ISSData.azimuth = Math.toDegrees(station.getAzimuth ...
    (ISSPosition, frame, currentDate));

idxs = ISSData.azimuth >= 180;
ISSData.azimuth(idxs) = ISSData.azimuth(idxs) - 360;

ISSData.elevation = Math.toDegrees(station.getElevation ...
    (ISSPosition, frame, currentDate));
ISSData.range = station.getRange(ISSPosition,...
    frame, currentDate);

end

%% Switch plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function id = getId(handles)
timeSeries = get(handles.TP, 'UserData');
t = timeSeries.t;

orekitData = get(handles.ISSTracker, 'UserData');
utc = orekitData.utc;

time = getCurrentTime(utc);
elapsed = time.durationFrom(orekitData.startTime);

id = find(elapsed > t, 1, 'last');
end

function switchPlots(handles)

set(handles.textS, 'Visible', 'off');
set(handles.textS2, 'Visible', 'off');
set(handles.textN, 'Visible', 'off');
set(handles.textE, 'Visible', 'off');
set(handles.textW, 'Visible', 'off');


timeSeries = get(handles.TP, 'UserData');
axes(handles.TP);
id = getId(handles);

t = timeSeries.t;

% newTicks = [fliplr(0:-3600:t(1)), 3600:3600:t(end)];

switch get(handles.whatplot, 'UserData')
    case 'elt'
        
        el = timeSeries.el;
        
        xlabel('Time (in s)');
        ylabel('Elevation (in °)');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', el);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', el(id));
        
        axis auto
        
    case 'azt'
        
        az = timeSeries.az;
        
        xlabel('Time (in s)');
        ylabel('Azimuth (in °)');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', az);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', az(id));
        
        axis auto
        
    case 'at'
        
        a = timeSeries.a;
        
        xlabel('Time (in s)');
        ylabel('Semi-major axis (in m)');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', a);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', a(id));
        
        axis auto
        
    case 'altt'
        
        alt = timeSeries.alt;
        
        xlabel('Time (in s)');
        ylabel('Altitude (in km)');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', alt);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', alt(id));
        
        axis auto
        
    case 'et'
        
        e = timeSeries.e;
        
        xlabel('Time (in s)');
        ylabel('Eccentricity');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', e);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', e(id));
        
        axis auto
        
    case 'it'
        
        i = timeSeries.i;
        
        xlabel('Time (in s)');
        ylabel('Inclination (in °)');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', i);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', i(id));
        
        axis auto
        
    case 'mat'
        
        ma = timeSeries.ma;
        
        xlabel('Time (in s)');
        ylabel('Mean anomaly (in °)');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', ma);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', ma(id));
        
        axis auto
        
    case 'pat'
        
        pa = timeSeries.pa;
        
        xlabel('Time (in s)');
        ylabel('Perigee Argument (in °)');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', pa);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', pa(id));
        
        axis auto
        
    case 'rat'
        
        raan = timeSeries.raan;
        
        xlabel('Time (in s)');
        ylabel('Right Ascension (in °)');
        
        set(handles.ISSTP(1), 'XData', t);
        set(handles.ISSTP(1), 'YData', raan);
        
        for k = 2:length(handles.ISSTP)
            set(handles.ISSTP(k), 'XData', []);
            set(handles.ISSTP(k), 'YData', []);
        end
        
        set(handles.ISSSta, 'XData', t(id));
        set(handles.ISSSta, 'YData', raan(id));
        
        axis auto
        
    case 'azel'
        az = timeSeries.az;
        el = timeSeries.el;
        
        myDif = abs(az(2:end) - az(1:end-1));
        idx = [0,find(myDif > 100),length(az)];
        
        for k=2:length(idx)
            set(handles.ISSTP(k-1), 'XData', az(idx(k-1)+1:idx(k)));
            set(handles.ISSTP(k-1), 'YData', el(idx(k-1)+1:idx(k)));
        end
        
        set(handles.ISSSta, 'XData', az(id));
        set(handles.ISSSta, 'YData', el(id));
        
        xlabel('Azimuth (in °)');
        ylabel('Elevation (in °)');
        
        axis([-180 180 -90 90])
        
        set(handles.textS, 'Visible', 'on');
        set(handles.textS2, 'Visible', 'on');
        set(handles.textN, 'Visible', 'on');
        set(handles.textE, 'Visible', 'on');
        set(handles.textW, 'Visible', 'on');
end

end

function switch3DPlots(handles)

hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

plot3D = get(handles.uipanel1, 'UserData');

switch plot3D
    case 'earth'
        set(handles.uipanel1, 'UserData', 'sky');
        setVisibilityThDV(handles, 'off');
        setVisibilitySkyMap(handles, 'on')
    case 'sky'
        set(handles.uipanel1, 'UserData', 'earth');
        setVisibilityThDV(handles, 'on');
        setVisibilitySkyMap(handles, 'off')
end


end

function setVisibilityThDV(handles, visible)
arrayfun(@(x)set(x, 'Visible', visible), handles.ThDVChildren);
end

function setVisibilitySkyMap(handles, visible)
arrayfun(@(x)set(x, 'Visible', visible), handles.skyMapChildren);
end


%% Used methods created by GUIDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes when startTrackerButton is pressed
function startTrackerButton_Callback(~, ~, ~)
% hObject    handle to startTrackerButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

value = get(handles.startTrackerButton,'UserData');
if ~value
    set(handles.startTrackerButton, 'String', 'Stop');
    set(handles.startTrackerButton, 'UserData', 1); % running
    drawnow;
    startTracking(hObject, handles);
else
    set(handles.startTrackerButton, 'String', 'Start');
    set(handles.startTrackerButton, 'UserData', 0); % not running
end

end

% --- Executes when user attempts to close ISSTracker.
function ISSTracker_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to ISSTracker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
set(handles.startTrackerButton, 'UserData', 0);
drawnow;
pause(.05);

infoFig = findall(0, 'Name', 'ISSTrackerInfo');
if ~isempty(infoFig)
    delete(infoFig)
end
pause(.2);

delete(hObject);

end

% --- Executes on button press in paltt.
function paltt_Callback(~, ~, handles)
% hObject    handle to paltt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'altt');
switchPlots(handles);
end

% --- Executes on button press in pazel.
function pazel_Callback(~, ~, handles)
% hObject    handle to pazel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'azel');
switchPlots(handles);
end

% --- Executes on button press in pmat.
function pmat_Callback(~, ~, handles)
% hObject    handle to pmat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'mat');
switchPlots(handles);
end

% --- Executes on button press in prat.
function prat_Callback(~, ~, handles)
% hObject    handle to prat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'rat');
switchPlots(handles);
end

% --- Executes on button press in ppat.
function ppat_Callback(~, ~, handles)
% hObject    handle to ppat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'pat');
switchPlots(handles);
end

% --- Executes on button press in pit.
function pit_Callback(~, ~, handles)
% hObject    handle to pit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'it');
switchPlots(handles);

end

% --- Executes on button press in pet.
function pet_Callback(~, ~, handles)
% hObject    handle to pet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'et');
switchPlots(handles);
end

% --- Executes on button press in pat.
function pat_Callback(~, ~, handles)
% hObject    handle to pat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'at');
switchPlots(handles);
end


% --- Executes on button press in pelt.
function pelt_Callback(~, ~, handles)
% hObject    handle to pelt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'elt');
switchPlots(handles);
end

% --- Executes on button press in pazt.
function pazt_Callback(~, ~, handles)
% hObject    handle to pazt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.whatplot, 'UserData', 'azt');
switchPlots(handles);
end

% --- Executes on button press in toggleSkyMap.
function toggleSkyMap_Callback(~, ~, ~)
% hObject    handle to toggleSkyMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch3DPlots;
% Hint: get(hObject,'Value') returns toggle state of toggleSkyMap
end

%% Unused methods created by GUIDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function timeInfo_Callback(hObject, eventdata, handles)
% hObject    handle to timeInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeInfo as text
%        str2double(get(hObject,'String')) returns contents of timeInfo as a double

end
function timeInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
function a_Callback(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a as text
%        str2double(get(hObject,'String')) returns contents of a as a double
end
function a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function e_Callback(hObject, eventdata, handles)
% hObject    handle to e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e as text
%        str2double(get(hObject,'String')) returns contents of e as a double
end
function e_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function i_Callback(hObject, eventdata, handles)
% hObject    handle to i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i as text
%        str2double(get(hObject,'String')) returns contents of i as a double
end
function i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
function pa_Callback(hObject, eventdata, handles)
% hObject    handle to pa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pa as text
%        str2double(get(hObject,'String')) returns contents of pa as a double
end
function pa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function raan_Callback(hObject, eventdata, handles)
% hObject    handle to raan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of raan as text
%        str2double(get(hObject,'String')) returns contents of raan as a double
end
function raan_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
% hObject    handle to raan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function ma_Callback(hObject, eventdata, handles)
% hObject    handle to ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ma as text
%        str2double(get(hObject,'String')) returns contents of ma as a double
end
function ma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function lat_Callback(hObject, eventdata, handles)
% hObject    handle to lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lat as text
%        str2double(get(hObject,'String')) returns contents of lat as a double
end
function lat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function lon_Callback(hObject, eventdata, handles)
% hObject    handle to lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lon as text
%        str2double(get(hObject,'String')) returns contents of lon as a double
end
function lon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function alt_Callback(hObject, eventdata, handles)
% hObject    handle to alt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alt as text
%        str2double(get(hObject,'String')) returns contents of alt as a double
end
function alt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function el_Callback(hObject, eventdata, handles)
% hObject    handle to el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of el as text
%        str2double(get(hObject,'String')) returns contents of el as a double
end
function el_CreateFcn(hObject, eventdata, handles)
% hObject    handle to el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function az_Callback(hObject, eventdata, handles)
% hObject    handle to az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of az as text
%        str2double(get(hObject,'String')) returns contents of az as a double
end
function az_CreateFcn(hObject, eventdata, handles)
% hObject    handle to az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function ra_Callback(hObject, eventdata, handles)
% hObject    handle to ra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ra as text
%        str2double(get(hObject,'String')) returns contents of ra as a double
end
function ra_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function vel_Callback(hObject, eventdata, handles)
% hObject    handle to vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vel as text
%        str2double(get(hObject,'String')) returns contents of vel as a double
end
function vel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function boundRotation_Callback(hObject, eventdata, handles)
% hObject    handle to boundRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boundRotation
end

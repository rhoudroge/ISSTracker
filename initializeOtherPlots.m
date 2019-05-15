% --- Initialization of telescope pointing axes
% Store other related time plots in UserData
function initializeOtherPlots(az, el, t, a, e, i, pa, ma, raan, alt, id)

hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

lat = simulation_parameters.userLat;
lon = simulation_parameters.userLon;
loc = simulation_parameters.userLoc;

str = 'Telescope Pointing - User location is';
if ~isempty(hObject)
    set(handles.uipanel4, 'Title', sprintf('%s %s (%.6f°, %.6f°)',...
        str, loc, lat, lon));
end

timeSeries.t = t;
timeSeries.az = az;
timeSeries.el = el;
timeSeries.a = a;
timeSeries.e = e;
timeSeries.i = i;
timeSeries.pa = pa;
timeSeries.ma = ma;
timeSeries.raan = raan;
timeSeries.alt = alt;

set(handles.TP, 'UserData', timeSeries);

axes(handles.TP);
hold on

% darkGray = [1 1 1] * .15;
xlabel('Azimuth (in °)');
ylabel('Elevation (in °)');

grid on

set(handles.TP, 'YColor', [.2 .2 .2]);
set(handles.TP, 'XColor', [.2 .2 .2]);

timeSeries = get(handles.TP, 'UserData');
az = timeSeries.az;
el = timeSeries.el;

myDif = abs(az(2:end) - az(1:end-1));
idx = [0,find(myDif > 100),length(az)];

for k=2:length(idx)
    handles.ISSTP(k-1) = plot(az(idx(k-1)+1:idx(k)), el(idx(k-1)+1:idx(k)), 'r');
end

handles.ISSSta = plot(az(id), el(id), 'y+', 'LineWidth', 2);
set(handles.TP, 'Color', 'k');

axis([-180 180 -90 90]);

guidata(hObject, handles);


end
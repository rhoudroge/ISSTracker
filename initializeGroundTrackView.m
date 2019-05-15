% --- Initialization of 2D View / Executes on load
%
% PARAMETERS
%   lons   array of longitudes
%   lats   array of latitudes
%   id     index of current date
%   r      radius at current date
function initializeGroundTrackView(lons, lats, id, r)

hObject = findall(0, 'Tag', 'ISSTracker');
handles = guidata(hObject);

axes(handles.ToDV); %#ok<*MAXES>

% Draw Earth and ISS on map
im = imread(fullfile('data', 'earthMapBM.jpg'));
image(0:360,-90:90,im);

axis equal
axis off

hold on

% adjust longitudes and latitudes
lons = lons + 180;
lats = -lats;

myDif = abs(lons(2:end) - lons(1:end-1));
idx = [0,find(myDif > 250),length(lons)];

% draw ground trace
%  avoid a line from 180° to -180°
for k=2:length(idx)
    plot(lons(idx(k-1)+1:idx(k)), lats(idx(k-1)+1:idx(k)), 'y','LineSmoothing','on');
end

%  arrow on actual ISS position
handles.myISS2D = plot(lons(id), lats(id), 'r+','LineWidth', 2,'LineSmoothing','on');

% footprint - r already normalized
beta = acos(1 / r);

handles.footprintArray = linspace(0, 2*pi, simulation_parameters.footPrint);
footLons = lons(id) + rad2deg(beta * cos(handles.footprintArray));
footLats = lats(id) + rad2deg(beta * sin(handles.footprintArray));

idxs = footLons > 360 | footLons < 0;

ffootLons = [];
ffootLats = [];

if sum(idxs) > 0
    
    if sum(footLons > 360) >= 1
        % right side
        ffootLons = footLons(idxs) - 360;
        n = length(ffootLons);
        ffootLons = [ffootLons(1:n/2) 0 0 ffootLons(n/2+1:n)];
        ll = rad2deg(beta) * sin(acos((360 - lons(id))/rad2deg(beta)));
        ll1 = -ll + lats(id);
        ll2 =  ll + lats(id);
        footLons(idxs) = 360;
        ffootLats = footLats(idxs);
        ffootLats = [ffootLats(1:n/2) ll2 ll1 ffootLats(n/2+1:n)];
    else
        % left side
        ffootLons = [360 footLons(idxs) + 360 360];
        ll = rad2deg(beta) * sin(acos(lons(id)/rad2deg(beta)));
        ll1 = -ll + lats(id);
        ll2 =  ll + lats(id);
        footLons(idxs) = 0;
        ffootLats = [ll2 footLats(idxs) ll1];
    end
    
end



handles.ISSFootPrint = fill(footLons, footLats, 'g', 'EdgeColor', 'g', ...
    'EdgeAlpha', .45, 'FaceColor', 'g', 'FaceAlpha', .2);
if isempty(ffootLons)
    handles.ISSFootPrintMinor = fill(0, 0, 'g', 'EdgeColor', ...
        'g', 'EdgeAlpha', .45, 'FaceColor', 'g', 'FaceAlpha', .2, 'Visible', 'off');
else
    handles.ISSFootPrintMinor = fill(ffootLons, ffootLats, 'g', 'EdgeColor', ...
        'g', 'EdgeAlpha', .45, 'FaceColor', 'g', 'FaceAlpha', .2, 'Visible', 'on');
end

guidata(hObject, handles);

end
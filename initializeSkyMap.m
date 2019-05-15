
% --- Initialization of 2D View / Executes on load
%
function initializeSkyMap(visible)

import org.orekit.time.*;

% locate ISSTracker
hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

% read catalog
fid=fopen(fullfile('data', simulation_parameters.skyMapCatalog));
M=textscan(fid, simulation_parameters.catalogRegEx, 'headerlines', 1);

% parse and keep stars over min magnitude
ct = getCurrentTime(TimeScalesFactory.getUTC);
c = ct.getComponents(TimeScalesFactory.getUTC).getDate.getYear;
YYYY = c(1);

sao=M{1};
mvis=M{2};
RA=(M{3}+(M{5}*(YYYY-2000)))/15;
DEC=(M{4}+(M{6}*(YYYY-2000)));

minMag = simulation_parameters.minMag;
inf4idxs = mvis < minMag;

sao = sao(inf4idxs);
mvis = mvis(inf4idxs);
RA = RA(inf4idxs);
DEC = DEC(inf4idxs);

starData{1} = sao;
starData{2} = RA;
starData{3} = DEC;

set(handles.skyMap, 'UserData', starData);

% select axes
axes(handles.skyMap);

% start drawing
hold on

blueGray = [103 219 244]/255;
darkBlueGray = [5 56 67]/255;

rectangle('Position', [-90 -90 180 180], 'Curvature', [1 1], ...
    'EdgeColor', darkBlueGray, 'FaceColor', 'k');

handles.skyLines(1) = plot([0 0 0],[-90 90 -90],'Color',darkBlueGray,'linewidth',1);
handles.skyLines(2) = plot([-90 90 -90],[0 0 0],'Color',darkBlueGray,'linewidth',1);


rs = linspace(30, 70, 3);
angs = linspace(0, 2*pi, 100);
ca = cos(angs);
sa = sin(angs);

for i = 1:length(rs)
    r = rs(i);
    xs = r * ca;
    ys = r * sa;
    plot(xs, ys,'Color',darkBlueGray,'LineSmoothing','on')
end

angs = [30 60 120 150 210 240 300 330];

rmi = min(rs);
rma = 90;

for i = 1:length(angs)
    ca = cos(angs(i) * pi / 180);
    sa = sin(angs(i) * pi / 180);
    xs = [rmi * ca, rma * ca];
    ys = [rmi * sa, rma * sa];
    plot(xs, ys,'Color',darkBlueGray,'LineSmoothing','on');
end

handles.star = zeros(length(sao), 1);
for i=1:length(sao)
    handles.star(i) = scatter(0,0,((minMag-mvis(i))*2),'*','filled','MarkerEdgeColor','white', 'MarkerFaceColor', 'white');
end

% iss trajectory
handles.ISSTrail1 = plot(0, 0, 'Color', blueGray, 'LineSmoothing','on',...
    'LineWidth', 1.5, 'Visible', 'off');
handles.ISSTrail2 = plot(0, 0, 'Color', blueGray, 'LineSmoothing','on',...
    'LineWidth', 1.5, 'Visible', 'off');
handles.ISSTrail3 = plot(0, 0, 'Color', blueGray, 'LineSmoothing','on',...
    'LineWidth', 1.5, 'Visible', 'off');
handles.ISSTrail4 = plot(0, 0, 'Color', blueGray, 'LineSmoothing','on',...
    'LineWidth', 1.5, 'Visible', 'off');
handles.ISSTrail5 = plot(0, 0, 'Color', blueGray, 'LineSmoothing','on',...
    'LineWidth', 1.5, 'Visible', 'off');

handles.ISSSky = plot(0, 0, 'Color', blueGray, 'Marker', '^', ...
    'LineSmoothing', 'on', 'LineWidth', 1, 'Visible', 'off');

axis equal tight;
axis off;
axis([-90 90 -90 90])

% make a list of all elements belonging to this plot
handles.skyMapChildren = [allchild(handles.skyMap); handles.star];

% upload new handles
guidata(hObject, handles);

% initial position of stars
updateSkyMapStars(TimeScalesFactory.getUTC);

end


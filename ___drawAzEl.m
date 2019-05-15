% draw the azimuth elevation graph
function drawAzEl(id)

hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

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
    handles.ISSTP(k-1) = plot(az(idx(k-1)+1:idx(k)), el(idx(k-1)+1:idx(k)), 'r','LineSmoothing','on');
end

handles.ISSSta = plot(az(id), el(id), 'y+', 'LineWidth', 2,'LineSmoothing','on');
set(handles.TP, 'Color', 'k');

axis([-180 180 -90 90]);

guidata(hObject, handles);

end
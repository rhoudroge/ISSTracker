% --- Initialization of 3D View / Executes on load
function initializeEarthView(x, y, z, id, lon, lat, ERA, zcomp, visible)
% Draw Earth and ISS in ITRF Frame

hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

axes(handles.ThDV);
cla

hold on

[handles.ellipsoid.X,handles.ellipsoid.Y,handles.ellipsoid.Z] =...
    ellipsoid(0,0,0,1,1,(1 - physical_parameters.f),30);

if simulation_parameters.renderEarth
    % draw map
    handles.im = imread(fullfile('data','earthMapBMsmall.jpg'));
    for k=1:3
        handles.im(:,:,k) = flipud(handles.im(:,:,k));
    end
    
    % Get image width in pixels
    width = size(handles.im, 2);
    
    % determine the angle per pixel
    anglePerPx = 2 * pi / width;
    
    % determine how many pixels need to be moved
    pxs = floor(ERA / anglePerPx);
    
    % Earth rotation!!
    for k=1:3
        current = handles.im(:,:,k);
        if zcomp > 0
            newCdata(:,:,k) = [current(:, width - pxs:width), current(:, 1:width - pxs - 1)];
        else
            newCdata(:,:,k) = [current(:, pxs + 1:width), current(:, 1:pxs)];
            %         newCdata(:,:,k) = [current(:, width - pxs:width), current(:, 1:width - pxs - 1)];
        end
    end
    
    handles.myEarth = surf(handles.ellipsoid.X, handles.ellipsoid.Y,...
        handles.ellipsoid.Z, 'Facecolor','texturemap',...
        'CData', newCdata, 'Edgecolor','none',...
        'FaceLighting','phong','LineSmoothing','on');
    
    
else
    % draw coastline only
    p1 = mesh(handles.ellipsoid.X, handles.ellipsoid.Y,...
        handles.ellipsoid.Z,'LineSmoothing','on');
    set(p1,'tag','earth','facecolor',[0 0 1],'edgecolor',[.3 .3 1]);
    
    data = load(fullfile('data','coastline'));
    coast = data.coastline;
    handles.coastLineTheta = coast(:,1)*pi/180;
    handles.coastLinePhi = coast(:,2)*pi/180;
    
    theta = handles.coastLineTheta + zcomp * ERA;
    phi = handles.coastLinePhi;
    coastLine = [cos(theta).*cos(phi),...
        sin(theta).*cos(phi),...
        -sin(phi)];
    handles.myEarth = plot3(coastLine(:,1),coastLine(:,2),-coastLine(:,3));
    set(handles.myEarth,'color',[0 .9 0]);
    
    
end

axis equal
axis off

% Plot EME2000 frame axes
h1 = plot3([0 1.5],[0 0],[0 0],'r-','LineSmoothing','on');
h2 = plot3([0 0],[0 1.5],[0 0],'g-','LineSmoothing','on');
h3 = plot3([0 0],[0 0],[0 1.5],'b-','LineSmoothing','on');
set([h1 h2 h3],'linewidth',4);

s = .05;
plotX = [0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*s -.025 + x(id);
plotY = [0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*s -.025 + y(id);
plotZ = [0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*s -.025 + z(id);

set(handles.ThDV, 'View', [lon(id) + 90 + sign(zcomp) * ERA * 180 / pi, ...
    lat(id) * simulation_parameters.latCoef]);

for i=1:6
    handles.cube(i) = patch(plotX(:,i),plotY(:,i),plotZ(:,i),'w');
    set(handles.cube(i),'edgecolor','k')
end

handles.orbit = plot3(x, y, z, 'y', 'LineWidth',1,'LineSmoothing','on');

% TODO create 3D footprint
x = x(id);
y = y(id);
z = z(id);
r = sqrt(x*x + y*y + z*z);

xr = x / r;
yr = y / r;
zr = z / r;

beta = acos(1 / r);

% list all children of ThDV axes
handles.ThDVChildren = allchild(handles.ThDV);
guidata(hObject, handles);

arrayfun(@(x)set(x, 'Visible', visible), handles.ThDVChildren);


end
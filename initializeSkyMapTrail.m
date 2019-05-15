function initializeSkyMapTrail(azs, els, times)

hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

azs0 = azs(els > 0);
els0 = els(els > 0);
times0 = times(els > 0);

azs1 = azs0 - 270;

if (~isempty(els0))
    % adjust to correct quadrant
    [issXs, issYs] = adjustQuadrant(azs1, els0);
    
    % Look for breaks (signaling different passes)
    diff = abs(issXs(2:end) - issXs(1:end-1));
    idx = find(diff > 50);
    
    if length(idx) == 1
        issXs1 = issXs(1:idx(1));
        issYs1 = issYs(1:idx(1));
        ts1 = times0(1:idx(1));
        
        issXs2 = issXs(idx(1)+1:end);
        issYs2 = issYs(idx(1)+1:end);
        ts2 = times0(idx(1)+1:end);
        
        set(handles.ISSTrail1, 'XData', issXs1, 'Visible', 'on');
        set(handles.ISSTrail1, 'YData', issYs1);
        
        set(handles.ISSTrail2, 'XData', issXs2, 'Visible', 'on');
        set(handles.ISSTrail2, 'YData', issYs2);
        
        createTickMarks(issXs1, issYs1, ts1);
        createTickMarks(issXs2, issYs2, ts2);
        
    elseif length(idx) == 2
        issXs1 = issXs(1:idx(1));
        issYs1 = issYs(1:idx(1));
        ts1 = times0(1:idx(1));
        
        issXs2 = issXs(idx(1)+1:idx(2));
        issYs2 = issYs(idx(1)+1:idx(2));
        ts2 = times0(idx(1)+1:idx(2));
        
        issXs3 = issXs(idx(2)+1:end);
        issYs3 = issYs(idx(2)+1:end);
        ts3 = times0(idx(2)+1:end);
        
        set(handles.ISSTrail1, 'XData', issXs1, 'Visible', 'on');
        set(handles.ISSTrail1, 'YData', issYs1);
        
        set(handles.ISSTrail2, 'XData', issXs2, 'Visible', 'on');
        set(handles.ISSTrail2, 'YData', issYs2);
        
        set(handles.ISSTrail3, 'XData', issXs3, 'Visible', 'on');
        set(handles.ISSTrail3, 'YData', issYs3);
        
        createTickMarks(issXs1, issYs1, ts1);
        createTickMarks(issXs2, issYs2, ts2);
        createTickMarks(issXs3, issYs3, ts3);
    elseif length(idx) >= 3
        issXs1 = issXs(1:idx(1));
        issYs1 = issYs(1:idx(1));
        ts1 = times0(1:idx(1));
        
        issXs2 = issXs(idx(1)+1:idx(2));
        issYs2 = issYs(idx(1)+1:idx(2));
        ts2 = times0(idx(1)+1:idx(2));
        
        issXs3 = issXs(idx(2)+1:idx(3));
        issYs3 = issYs(idx(2)+1:idx(3));
        ts3 = times0(idx(2)+1:idx(3));
        
        issXs4 = issXs(idx(3)+1:end);
        issYs4 = issYs(idx(3)+1:end);
        ts4 = times0(idx(3)+1:end);
        
        set(handles.ISSTrail1, 'XData', issXs1, 'Visible', 'on');
        set(handles.ISSTrail1, 'YData', issYs1);
        
        set(handles.ISSTrail2, 'XData', issXs2, 'Visible', 'on');
        set(handles.ISSTrail2, 'YData', issYs2);
        
        set(handles.ISSTrail3, 'XData', issXs3, 'Visible', 'on');
        set(handles.ISSTrail3, 'YData', issYs3);
        
        set(handles.ISSTrail4, 'XData', issXs4, 'Visible', 'on');
        set(handles.ISSTrail4, 'YData', issYs4);
        
        createTickMarks(issXs1, issYs1, ts1);
        createTickMarks(issXs2, issYs2, ts2);
        createTickMarks(issXs3, issYs3, ts3);
        createTickMarks(issXs4, issYs4, ts4);
    else
        set(handles.ISSTrail1, 'XData', issXs, 'Visible', 'on');
        set(handles.ISSTrail1, 'YData', issYs);
        
        createTickMarks(issXs, issYs, times0);
    end
end

end

function [issXs, issYs] = adjustQuadrant(azs, els)

issXs = zeros(size(azs));
issYs = zeros(size(azs));

caz = cosd(azs);
saz = sind(azs);

% Adjust trail to correct quadrant
for i = 1:length(caz)
    
    ca = caz(i);
    sa = saz(i);
    
    [xxx,yyy]=linecircme(sa/ca,0,0,0,90-els(i));
    
    if (sa>=0) && (ca<=0)
        %1st Quadrant anticlockwise from North
        issXs(i) = min(xxx);
        issYs(i) = max(yyy);
    elseif (ca<=0) && sa<0
        %2nd Quadrant anticlockwise from North
        issXs(i) = min(xxx);
        issYs(i) = min(yyy);
    elseif ca>0 && (sa<=0)
        %3rd Quadrant anticlockwise from North
        issXs(i) = max(xxx);
        issYs(i) = min(yyy);
    elseif ca>0 && sa>0
        %4th Quadrant anticlockwise from North
        issXs(i) = max(xxx);
        issYs(i) = max(yyy);
    end
end

end


function createTickMarks(x, y, times)
import org.orekit.time.*;

% nombre et taille des marqueurs
n = simulation_parameters.tickMarks;
l = 2;

% tableau de secondes écoulées depuis le début de times
t = cellfun(@(xi) xi.durationFrom(times{1}), times);

% coordonnées des marqueurs
mx = linspace(x(1), x(end), n+2);
mx = mx(2:end-1);
my = arrayfun(@(xi) interp1(x, y, xi), mx);
mt = arrayfun(@(xi) interp1(x, t, xi), mx);
dates = cell(length(mt),1);
utc = TimeScalesFactory.getUTC;
for i=1:length(mt)
    dates{i} = char(times{1}.shiftedBy(mt(i)).toString(utc));
end

% abscisses + et - pour calcul de la tangeante
mxp = mx + l;
mxm = mx - l;

% ordonnées + et - pour calcul de la tangeante
myp = arrayfun(@(xi) interp1(x, y, xi), mxp);
mym = arrayfun(@(xi) interp1(x, y, xi), mxm);

% coefficients de la tangeante
a_p = arrayfun(@(xp, xm, yp, ym) (yp - ym)/(xp - xm), mxp, mxm, myp, mym);

% coefficients de la perpendiculaire à la tangeante
a1_p = -1 ./ a_p;
b1_p = my - a1_p .* mx;

% abscisses des points de la perpendiculaire
mxpp = mx + sqrt(l ./ (1 + a1_p .* a1_p));
mxpm = mx - sqrt(l ./ (1 + a1_p .* a1_p));

mp_p = a1_p .* mxpp + b1_p;
mp_m = a1_p .* mxpm + b1_p;

% draw figure
hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);
axes(handles.skyMap);
blueGray = [103 219 244]/255;

for i=1:n
    plot([mxm(i) mxp(i)], [mp_m(i) mp_p(i)], 'Color', blueGray, ...
        'LineWidth', 1);
    text(mx(i), my(i)+5, dates{i}, 'Color', blueGray, 'FontSize', 7);
end

% locate ISSTracker
hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

% make a list of all elements belonging to this plot
handles.skyMapChildren = allchild(handles.skyMap);

% upload new handles
guidata(hObject, handles);

end

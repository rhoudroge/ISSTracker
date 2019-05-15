function updateSkyMapISS(az, el)

hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);

if (el > 0)
    ca = cosd(az - 270);
    sa = sind(az - 270);
    
    [xxx,yyy]=linecirc(sa/ca,0,0,0,90-el);
    
    if (sa>=0) && (ca<=0)
        %1st Quadrant anticlockwise from North
        issXs = min(xxx);
        issYs = max(yyy);
    elseif (ca<=0) && sa<0
        %2nd Quadrant anticlockwise from North
        issXs = min(xxx);
        issYs = min(yyy);
    elseif ca>0 && (sa<=0)
        %3rd Quadrant anticlockwise from North
        issXs = max(xxx);
        issYs = min(yyy);
    elseif ca>0 && sa>0
        %4th Quadrant anticlockwise from North
        issXs = max(xxx);
        issYs = max(yyy);
    end
    
    set(handles.ISSSky, 'Visible', 'on');
    set(handles.ISSSky, 'XData', issXs);
    set(handles.ISSSky, 'YData', issYs);
else
    set(handles.ISSSky, 'XData', 0);
    set(handles.ISSSky, 'YData', 0);
    set(handles.ISSSky, 'Visible', 'off');
end

end
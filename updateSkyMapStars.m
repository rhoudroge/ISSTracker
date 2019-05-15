% --- Update star positions on sky map
function updateSkyMapStars(utc)

import org.orekit.time.*;

currentTime = getCurrentTime(utc);

comps = currentTime.getComponents(utc);
ddate = comps.getDate;
ttime = comps.getTime;

YYYY = ddate.getYear;
MM   = ddate.getMonth;
DD   = ddate.getDay;
hh   = ttime.getHour;
mm   = ttime.getMinute;
ss   = ttime.getSecond;

% user location - degrees
lat = simulation_parameters.userLat;
long = simulation_parameters.userLon;

% star data
hObject = findall(0,'Tag','ISSTracker');
handles = guidata(hObject);
starData = get(handles.skyMap, 'UserData');
sao = starData{1};
RA = starData{2};
DEC = starData{3};

%-----------------------------------------------------------------
% GREENWICH MEAN SIDEREAL TIME
%-----------------------------------------------------------------
dfrac=(hh+mm/60+ss/3600)/24;
dwhole =367*YYYY-fix(7*(YYYY+fix((MM+9)/12))/4)+fix(275*MM/9)+DD-730531.5;
dd=dwhole+dfrac;  %no. of days elapsed since J2000.0
% dd=d-UTdiff/24;  %no of UT days elapsed since UT J2000.0
T=dd/36525;  %fraction of epoch/century time elapsed

% Meeus formula 11.4 for mean sidereal time at zero longitude (Greenwich Mean Sidereal Time)
GMST=280.46061837+(360.98564736629*dd)+(0.000387933*T^2)-(T^3/38710000);
while GMST>360
    GMST=GMST-360;
end
while GMST<0
    GMST=GMST+360;
end

%------------------------------------------------------------------
% LOCAL SIDEREAL TIME
%------------------------------------------------------------------
LST=GMST+long;
while LST>360
    LST=LST-360;
end
while LST<0
    LST=LST+360;
end
LST=LST/15;  %LST in hours

for i=1:length(sao)
    
    HA =(LST - RA(i))*15;
    El = asind(sind(lat)*sind(DEC(i))+cosd(lat)*cosd(DEC(i))*cosd(HA));  %http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
    
    if El>0   %star is in horizon
        if El~=90 %star not at zenith
            
            AZ1 = (-cosd(DEC(i))*sind(HA))/cosd(El);
            AZ2 = ( cosd(lat) * sind(DEC(i)) - sind(lat) * cosd(DEC(i)) * cosd(HA) ) / cosd(El);
            
            [xxx,yyy]=linecircme(AZ2/AZ1,0,0,0,90-El);
            
            if (AZ1>=0) && (AZ2>=0)
                %1st Quadrant anticlockwise from North
                set(handles.star(i), 'XData', min(xxx));
                set(handles.star(i), 'YData', max(yyy));
            elseif (AZ1>=0) && AZ2<0
                %2nd Quadrant anticlockwise from North
                set(handles.star(i), 'XData', min(xxx));
                set(handles.star(i), 'YData', min(yyy));
            elseif AZ1<0 && (AZ2>=0)
                %4th Quadrant anticlockwise from North
                set(handles.star(i), 'XData', max(xxx));
                set(handles.star(i), 'YData', max(yyy));
            elseif AZ1<0 && AZ2<0
                %3rd Quadrant anticlockwise from North
                set(handles.star(i), 'XData', max(xxx));
                set(handles.star(i), 'YData', min(yyy));
                
            end
        elseif El==90 %star at zenith
            %                 las = scatter(0,0,((maxmag-mvis(i))*2),'*','filled','MarkerEdgeColor','white', 'MarkerFaceColor', 'white'); hold on;
        end
    else
        set(handles.star(i), 'Visible', 'off');
    end
end


end

function skymapf

clear all
close all
clc

% YYYY=2014;
% MM=2;
% DD=16;
% hh=13;
% mm=0;
% ss=0;
% long=7.3;
% lat=6.8;
% UTdiff=+1;

import java.lang.System;
import java.io.File;
import org.orekit.data.*;
import org.orekit.time.*;

% data pointer
DM=DataProvidersManager.getInstance();
crawler=DirectoryCrawler(File(fullfile(cd, 'data')));
DM.clearProviders();
DM.addProvider(crawler);

maxmag = 5;

% timeNow = getTimeComponents;
utc = TimeScalesFactory.getUTC;
timeNow = getCurrentTime(utc);

YYYY = timeNow.getComponents(utc).getDate.getYear;
MM   = timeNow.getComponents(utc).getDate.getMonth;
DD   = timeNow.getComponents(utc).getDate.getDay;
hh   = timeNow.getComponents(utc).getTime.getHour;
mm   = timeNow.getComponents(utc).getTime.getMinute;
ss   = timeNow.getComponents(utc).getTime.getSecond;

% user location - degrees
lat = simulation_parameters.userLat;
long = simulation_parameters.userLon;

fid=fopen(fullfile('data', 'sao6.csv'));

M=textscan(fid, '%d;%f;%f;%f;%f;%f', 'headerlines', 1);

sao=M{1};
mvis=M{2};
RA=(M{3}+(M{5}*(YYYY-2000)))/15;
DEC=(M{4}+(M{6}*(YYYY-2000)));

inf4idxs = mvis < 4;

length(sao)
sao = sao(inf4idxs);
mvis = mvis(inf4idxs);
RA = RA(inf4idxs);
DEC = DEC(inf4idxs);
length(sao)


scrnsz=get(0,'ScreenSize');

pl = figure('Position',[scrnsz(1)+scrnsz(3)*0.1 scrnsz(2)+scrnsz(4)*0.1 scrnsz(3)*0.8 scrnsz(4)*0.8]);
set(pl, 'Color', [.7 .7 .7]);

hold on



rect = rectangle('Position', [-90 -90 180 180], 'Curvature', [1 1], ...
    'EdgeColor', 'g', 'FaceColor', 'k');

star = zeros(length(sao), 1);

for i=1:length(sao)
    star(i) = scatter(0,0,((maxmag-mvis(i))*2),'*','filled','MarkerEdgeColor','white', 'MarkerFaceColor', 'white');
    set(star(i), 'Visible', 'off');
end

title(['Sky Map for Longitude ' num2str(long) ', Latitude ' num2str(lat) ', ' num2str(YYYY) '-' num2str(MM) '-' num2str(DD) ', ' num2str(hh) ':' num2str(mm) ':' num2str(ss) ' UTC']);
ylabel(''); xlabel('');
axis equal tight;
axis off;
axis([-90 90 -90 90])





%% Time vector
import org.orekit.time.*;
import org.orekit.time.*;

utc = TimeScalesFactory.getUTC;
currentTime =  AbsoluteDate(DateComponents(YYYY, MM, DD),...
    TimeComponents(hh, mm, ss), utc);


%% update

seconds = 0:150:3600*4;

for iii = 1:length(seconds)
    
    myTime = currentTime.shiftedBy(seconds(iii));
    comps = myTime.getComponents(utc);
    ddate = comps.getDate;
    ttime = comps.getTime;
    
    YYYY = ddate.getYear;
    MM   = ddate.getMonth;
    DD   = ddate.getDay;
    hh   = ttime.getHour;
    mm   = ttime.getMinute
    ss   = ttime.getSecond
    
    
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
            set(star(i), 'Visible', 'on');
            if El~=90 %star not at zenith
                
                AZ1 = (-cosd(DEC(i))*sind(HA))/cosd(El);
                AZ2 = (cosd(lat)*sind(DEC(i))-sind(lat)*cosd(DEC(i))*cosd(HA))/cosd(El);
                
                [xxx,yyy]=linecirc(AZ2/AZ1,0,0,0,90-El);
                
                if (AZ1>=0) && (AZ2>=0)
                    %1st Quadrant anticlockwise from North
                    set(star(i), 'XData', min(xxx));
                    set(star(i), 'YData', max(yyy));
                elseif (AZ1>=0) && AZ2<0
                    %2nd Quadrant anticlockwise from North
                    set(star(i), 'XData', min(xxx));
                    set(star(i), 'YData', min(yyy));
                elseif AZ1<0 && (AZ2>=0)
                    %4th Quadrant anticlockwise from North
                    set(star(i), 'XData', max(xxx));
                    set(star(i), 'YData', max(yyy));
                elseif AZ1<0 && AZ2<0
                    %3rd Quadrant anticlockwise from North
                    set(star(i), 'XData', max(xxx));
                    set(star(i), 'YData', min(yyy));
                end
            end
        else
            set(star(i), 'Visible', 'off');
        end
    end
    
    title(['Sky Map for Longitude ' num2str(long) ', Latitude ' num2str(lat) ', ' num2str(YYYY) '-' num2str(MM) '-' num2str(DD) ', ' num2str(hh) ':' num2str(mm) ':' num2str(ss) ' LT']);
    drawnow
    
    
    
end

end


function timeNow = getTimeComponents

import java.lang.Math;
import java.lang.System;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.TimeZone;

% - Date of now
date = Date();
calendar = GregorianCalendar.getInstance(TimeZone.getTimeZone('UTC')); % creates a new calendar instance
calendar.setTime(date);   % assigns calendar to given date

year = calendar.get(Calendar.YEAR);
month = calendar.get(Calendar.MONTH)+1;
day = calendar.get(Calendar.DAY_OF_MONTH);
hour = calendar.get(Calendar.HOUR_OF_DAY); % gets hour in 24h format
minute = calendar.get(Calendar.MINUTE);
second = calendar.get(Calendar.SECOND) + calendar.get(Calendar.MILLISECOND)/1000 ;

timeNow = [year, month, day, hour, minute, second];
end

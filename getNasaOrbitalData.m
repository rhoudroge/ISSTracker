function orbitalData = getNasaOrbitalData( data )
% GETNASAORBITALPOSITIONS Returns the orbital data from the given data feed
%   orbitalData = getNasaOrbitalData(data) returns an orbitalData structure
%   based on the given data.
% See also getNasaActualFeed getNasaSavedFeed

% useful for DateComponents
import org.orekit.time.*;

% conversion from LBS to KGS
lbs2kgs = 0.453592;

% Split the data char starting from the startStr
startStr = '    Coasting Arc';
startIdx = [strfind(data, startStr) length(data)];

% For each state vector from NASA, retrieve the info we want
for k = 1:length(startIdx)-1
    % parse into lines
    orbitalData(k).text.raw = data(startIdx(k):startIdx(k+1)-1);
    orbitalData(k).text.lines = regexp(orbitalData(k).text.raw, '\n', 'split');
    
    % get time info
    timeData = orbitalData(k).text.lines{4};
    dateData = char(regexpi(timeData,...
        '\d+/\d+/\d\d:\d\d:\d\d.\d\d\d','match'));
    
    orbitalData(k).time.year = str2double(dateData(1:4));
    
    date = DateComponents(orbitalData(k).time.year, str2double(dateData([6 7 8])));
    
    orbitalData(k).time.month = date.getMonth;
    orbitalData(k).time.day = date.getDay;
    orbitalData(k).time.hour = str2double(dateData([10 11]));
    orbitalData(k).time.minute = str2double(dateData([13 14]));
    orbitalData(k).time.second = str2double(dateData([16 17]));
    orbitalData(k).time.msecond = str2double(dateData([19 20 21]));
    
    % get cartesian info
    labels = {'X', 'Y', 'Z', 'XDot', 'YDot', 'ZDot'};
    for j = 22:27
        result = regexp(orbitalData(k).text.lines{j},'(-?[0-9]{0,}.\d+)*','match');
        orbitalData(k).data.(labels{j-21}) = str2double(result{2});
    end
    
    orbitalData(k).data.weight = lbs2kgs * str2double( ...
        regexpi(orbitalData(k).text.lines{6}, '\d+.\d+', 'match'));
end


%     % Look for maneuvers
%     flagStr = '   IMPULSIVE TIG';
%     flag = strfind(nasaData, flagStr);
%     flagStr = '   IMPULSIVE TIG';
%     endFlag = strfind(nasaData, flagStr);
%     
%     maneuver = nasaData(flag(1):endFlag(1)-10);
%     
%     maneuverData.text.raw = maneuver;
%     maneuver = regexp(maneuver, '\n', 'split');
%     maneuverData.text.lines = maneuver;
%     dateData = char(regexpi(maneuver{5},...
%         '\d+/\d\d:\d\d:\d\d.\d\d\d','match'));
%     % if length(maneuver) > 6
%     %     maneuverData.start = AbsoluteDate(DateComponents(cal.get(Calendar.YEAR),...
%     %         str2double(dateData(1:3))), TimeComponents(str2double(dateData(5:6)), ...
%     %         str2double(dateData(8:9)), str2double(dateData(11:12)) + ...
%     %         str2double(dateData(14:16))/1000), utc);
%     %     dateData = char(regexpi(maneuver{7},...
%     %         '\d\d:\d\d:\d\d.\d\d\d           \d.\d','match'));
%     %     maneuverData.duration = str2double(dateData(1:2)) * 3600 + ...
%     %         str2double(dateData(4:5)) * 60 + str2double(dateData(7:8)) + ...
%     %         str2double(dateData(10:12)) / 1000;
%     %     maneuverData.containsManeuver = true;
%     % else
%     maneuverData.containsManeuver = false;
%     % end

end


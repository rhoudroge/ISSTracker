function saveNasaFeed
% saveNasaFeed Saves the current NASA ephemeris to the work space

import org.orekit.time.*;

% get the current time
configureLibraries;
utc = TimeScalesFactory.getUTC;
timeNow = getCurrentTime(utc, 1);
time = date2char(timeNow);

% create filename
filename = fullfile(cd, 'in', ['nasa_feed_' time(1:8) '.txt']);

% print to file
urlwrite(simulation_parameters.url, filename);

fprintf('Saved data to %s\n', filename);

end


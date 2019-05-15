function nasaData = getNasaActualFeed
% GETNASAACTUALFEED Returns the current NASA ephemeris as a Matlab String
%   See also getNasaOrbitalData.

url = simulation_parameters.url;
try
    nasaData = urlread(url);
catch e
    error('Couldn''t retrieve orbital data. Please check your internet connection.');
end

end


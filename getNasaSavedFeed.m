function nasaData = getNasaSavedFeed( filename )
% GETNASASAVEDFEED Returns the saved NASA file as a Matlab String
% nasaData = getNasaSavedFeed( filename ) returns the contents of filename.
% It is assumed that the given file is in the folder 'in'.
%   See also getNasaOrbitalData saveNasaFeed

nasaData = fileread(fullfile('in', filename));

end


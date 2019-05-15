function dateChar = date2char( date )
%DATETOCHAR Summary of this function goes here
%   Detailed explanation goes here

t0 = char(date.toString());
t0 = strrep(t0, ':', '');
t0 = strrep(t0, '-', '');
dateChar = t0(1:end-4);

end


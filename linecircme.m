function [xout, yout] = linecircme(a, b, x0, y0, r)
% [xout,yout] = linecirc(slope,intercpt,centerx,centery,radius) finds the
% points of intersection given a circle defined by a center and radius in
% x-y coordinates, and a line defined by slope and y-intercept, or a slope
% of "inf" and an x-intercept. Two points are returned. When the objects do
% not intersect, NaNs are returned.
%
% When the line is tangent to the circle, two identical points are
% returned. All inputs must be scalars.
%
% It does exactly the same thing as linecirc

% Given P(x, y), P belongs to the circle and the line if
%    y = a * x + b
%    A * x^2 + B * x + C = 0
A = 1 + a * a;
B = 2 * (a*(b - y0) - x0);
C = x0 * x0 + (b-y0)^2 - r * r;

% simplified delta = B^2 - 4 * A * C
% delta = 4 * (x02 - 2 * a * b * x0 + 2 * a * y0 * x0 ...
%     - x02 + r2 - b2 - a2 * x02 + a2 * r2 + 2 * b * y0 - y02);
delta = B * B - 4 * A * C;

if delta > 0
    
    xout(1) = (-B + sqrt(delta)) / (2 * A);
    yout(1) = a * xout(1) + b;
    
    xout(2) = (-B - sqrt(delta)) / (2 * A);
    yout(2) = a * xout(2) + b;
    
elseif delta == 0
    
    xout(1) = -B / (2 * A);
    yout(1) = a * xout(1) + b;
    
    xout(2) = xout(1);
    yout(2) = yout(1);
    
else
    xout = [NaN NaN];
    yout = xout;
    
end



end
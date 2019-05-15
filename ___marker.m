% fonction
x = 1:10;
y = x .* x / 10 + 1;

% nombre et taille des marqueurs
n = 7;
l = .2;

% coordonnées des marqueurs
mx = linspace(x(1), x(end), n);
my = arrayfun(@(xi) interp1(x, y, xi), mx);

% abscisses + et - pour calcul de la tangeante
mxp = mx + l;
mxm = mx - l;

% ordonnées + et - pour calcul de la tangeante
myp = arrayfun(@(xi) interp1(x, y, xi), mxp);
mym = arrayfun(@(xi) interp1(x, y, xi), mxm);

% coefficients de la tangeante
a_p = arrayfun(@(xp, xm, yp, ym) (yp - ym)/(xp - xm), mxp, mxm, myp, mym);
b_p = my - a_p .* mx;

% coefficients de la perpendiculaire à la tangeante
a1_p = -1 ./ a_p;
b1_p = my - a1_p .* mx;

% abscisses des points de la perpendiculaire
mxpp = mx + sqrt(l ./ (1 + a1_p .* a1_p));
mxpm = mx - sqrt(l ./ (1 + a1_p .* a1_p));

mp_p = a1_p .* mxpp + b1_p;
mp_m = a1_p .* mxpm + b1_p;

% draw figure
figure
plot(x, y, '-')
hold
for i=1:n
    plot([mxm(i) mxp(i)], [mp_m(i) mp_p(i)], 'r')
end
axis equal
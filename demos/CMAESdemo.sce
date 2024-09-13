// Author: Yann COLLETTE, 2006. 
// e-mail: yann[dot]colletet[at]renault[dot]com
// URL: http://ycollette.free.fr

// TestCMAES.sce

PlotMin   = %F;
PlotEigen = %F;

lines(0);

functoplot=rosenbrock;
Title='rosenbrock';
get_min=get_min_bound_rosenbrock;
get_max=get_max_bound_rosenbrock;
delta = 0.05;

x_min  = get_min();
x_max  = get_max();
x0     = (x_max - x_min).*rand(size(x_max,1),size(x_max,2)) + x_min;
x_step = delta*(x_max-x_min);

[xmin, fmin, _log] = cmaes(fsphere, x0, Log = %T);
PlotTitle = sprintf('%s function - xmin = [%f %f] - fmin = %f', Title, xmin(1), xmin(2), fmin);
plot_cmaes(_log,PlotTitle,plotmin=PlotMin, plot_eigen=PlotEigen);
[X, Y] = meshgrid(x_min(1):x_step(1):x_max(1), x_min(2):x_step(2):x_max(2));
for i=1:size(X,1)
  for j=1:size(X,2)
    Z(i,j) = functoplot([X(i,j);Y(i,j)]);
  end
end
scf(); surf(X,Y,Z);
xtitle(Title + ' Function','X','Y','Z');




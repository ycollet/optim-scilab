function [hs, x_history] = ls_polynom(polynom_f,polynom_fs,polynom_fss,x0,d,h,Log,TOL,ItMX,cond_func)
// polynomial method for f' 
//
// [hs, [x_history]]=ls_polynom(polynom_f,polynom_fs,polynom_fss,x0,h,Log,TOL,ItMX,cond_func)
//
// polynom_f    - handle to function
// polynom_fs   - handle to function evaluating f'
// polynom_fss  - handle to function evaluating f'' (ignored)
// x0   - initial guess
// h    - initial interval length
// Log  - print log information, 1=on (default), 0=off
// TOL  - convergence tolerance: (x3-x1)<TOL, default: 1.0e-8
// ItMX - maximal number of iterations, default: 25
// x_history - list of computed points in the parameter space

[nargout, nargin] = argn();

if (nargin < 6) then Log  = 1;                 end
if (nargin < 7) then TOL  = 1.0e-8;            end
if (nargin < 8) then ItMX = 25;                end

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end

y0  = f(x0);
ys0 = d'*fs(x0);

// First part: localization of the optimum

if (x_history_defined) then
  [h_list_aux, x_history_aux] = ls_2pts_bracket(polynom_f,polynom_fs,polynom_fss,x0,d,h,Log,TOL,ItMX);
  hk_min = h_list_aux(1);
  hk_max = h_list_aux(2);
  x_history = lstcat(x_history, x_history_aux);
else
  [h_list_aux] = ls_2pts_bracket(polynom_f,polynom_fs,polynom_fss,x0,d,h,Log,TOL,ItMX);
  hk_min = h_list_aux(1);
  hk_max = h_list_aux(2);
end

xk_min = x0 + hk_min*d;
xk_max = x0 + hk_max*d;

yk_min  = f(xk_min);
yk_max  = f(xk_max);
ysk_min = d'*polynom_fs(xk_min);

if (Log) then
  printf('ls_polynom\n');
  printf('%3s %10s %10s %10s %10s\n', 'It', '|d|', 'x', 'f(x)', 'fs(x)');
  printf('%3d %10.2e %10.2e %10.2e %10.2e\n', 0, 0.0, xk_min, yk_min, ysk_min);
end

// Second part: polynom identification

if (x_history_defined) then
  x_history($+1) = xk_min;
  x_history($+1) = xk_max;
end

a2 = ((yk_max - yk_min)/(hk_max - hk_min) - ysk_min)/(hk_max - hk_min);
a1 = ysk_min - 2*a2*hk_min;
a0 = yk_min - a1*hk_min - a2*hk_min^2;

b = a1^2 -4*a0*a2;
h1_sol = abs(real((-a1 + sqrt(b))/(2*a2)));
h2_sol = abs(real((-a1 - sqrt(b))/(2*a2)));

if (polynom_f(x0+h1_sol*d)<polynom_f(x0+h2_sol*d)) then
  hs = h1_sol;
else
  hs = h2_sol;
end

if (x_history_defined) then
  x_history($+1) = x0+hs*d;
end

if (Log) then
  printf('%3d %10.2e %10.2e %10.2e %10.2e\n', 1, hs, x0+hs*d, polynom_f(x0+hs*d), d'*polynom_fs(x0+hs*d));
end

// We transform hs as a list for coherency with other methods
aux     = hs;
hs      = list();
hs($+1) = aux;
endfunction

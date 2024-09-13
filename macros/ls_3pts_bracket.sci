function [h_list, x_history]=ls_3pts_bracket(ls_3pts_f,ls_3pts_fs,ls_3pts_fss,x0,d,h,Log,TOL,ItMX)
// method of 3 points bracketing for f
//
// [h_list,[x_history]]=ls_3pts_bracket(f,fs,fss,x0,h,Log,TOL,ItMX)
//
// ls_3pts_f   - handle to function
// ls_3pts_fs  - handle to function evaluating f'
// ls_3pts_fss - handle to function evaluating f'' (not used)
// x0   - initial guess
// h    - initial interval length
// Log  - print log information, 1=on (default), 0=off
// TOL  - convergence tolerance: (x3-x1)<TOL, default: 1.0e-8
// ItMX - maximal number of iterations, default: 25
// cond_func - function to test the validity of the step
// Output:
// h_history - list of computed points in the parameter space
// h_list - a list of step size: h_list(1) = h_min, h_list(2) = h_mid, h_list(3) = h_max

[nargout, nargin] = argn();

if (nargin < 6) then Log  = 1;                 end
if (nargin < 7) then TOL  = 1.0e-8;            end
if (nargin < 8) then ItMX = 25;                end

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end

h_min = 0;
h_mid = 0;
h_max = 0;

// Localization of the optimum

incr = 2;

y0  = ls_3pts_f(x0);
ys0 = d'*ls_3pts_fs(x0);

hk_max = h;
hk_min = 0.0;
hk_mid = (hk_max+hk_min)/2.0;

xk_max = x0 + d*hk_max;
xk_mid = x0 + d*hk_mid;
xk_min = x0 + d*hk_min;

yk_max = ls_3pts_f(xk_max);
yk_mid = ls_3pts_f(xk_mid);
yk_min = ls_3pts_f(xk_min);

if (x_history_defined) then
  x_history($+1) = xk_max;
  x_history($+1) = xk_mid;
  x_history($+1) = xk_min;
end

if (d'*ls_3pts_fs(xk_min)>0) then
  if Log then
    printf('ls_3pts_bracket: warning - initial direction descent not a direction descent\n');
  end
  h_mid = 0;
  h_min = 0;
  h_max = 0;
  
  h_list    = list();
  h_list(1) = h_min;
  h_list(2) = h_mid;
  h_list(3) = h_max;
  
  return;
end

if (Log) then
  printf('ls_3pts_bracket:\n');
  printf('%3s %10s %10s %10s %10s\n', 'It', 'h', 'f_min(h)', 'f_mid(h)', 'f_max(h)');
  printf('%3d %10.2e %10.2e %10.2e %10.2e\n', 0, x0, yk_min, yk_mid, yk_max);
end

// Looking for a good starting h_max: a h_max such that fs(x0+h_max*d)<0
It_start = 0;
ysk_max  = ls_3pts_fs(xk_max);
while((d'*ysk_max<0)&(It_start<ItMX))
  It_start = It_start + 1;
  if (x_history_defined) then
    x_history($+1) = xk_max;
  end

  hk_max = hk_max * incr;
  xk_max = x0 + d*hk_max;

  if (Log) then 
    yk_max = ls_3pts_f(xk_max);
    printf('%3d %10.2e %10.2e %10.2e %10.2e\n', It_start, hk_max, yk_min, yk_mid, yk_max);
  end
  
  ysk_max = ls_3pts_fs(xk_max);
end

hk_mid = (hk_max + hk_min) / 2;

h_min = hk_min;
h_mid = hk_mid;
h_max = hk_max;

h_list    = list();
h_list(1) = h_min;
h_list(2) = h_mid;
h_list(3) = h_max;

if (Log) then
  printf('ls_3pts_bracket: terminated after %d iterations\n', It_start);
end
endfunction

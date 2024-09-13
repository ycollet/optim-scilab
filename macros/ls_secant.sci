function [hs, x_history] = ls_secant(secant_f,secant_fs,secant_fss,x0,d,h,Log,TOL,ItMX)
// method of secants for f' [simplified 1d-Newton]
// not globalised
//
// [hs, [x_history]]=ls_secant(secant_f,x0,h,Log,TOL,ItMX)
//
// secant_f    - handle to function
// secant_fs   - handle to function evaluating f'
// secant_fss  - handle to function evaluating f'' (ignored)
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

// First part: localization of the optimum

if (x_history_defined) then
  [h_list_aux, x_history_aux] = ls_2pts_bracket(secant_f,secant_fs,secant_fss,x0,d,h,Log,TOL,ItMX);
  hk_min = h_list_aux(1);
  hk_max = h_list_aux(2);
  x_history = lstcat(x_history, x_history_aux);
else
  [h_list_aux] = ls_2pts_bracket(secant_f,secant_fs,secant_fss,x0,d,h,Log,TOL,ItMX);
  hk_min = h_list_aux(1);
  hk_max = h_list_aux(2);
end

// Second part: the secant process
y0  = secant_f(x0);
ys0 = secant_fs(x0);

xk    = x0+hk_min*d;
xk_1  = x0+hk_max*d;
xknew = xk;
ysk     = secant_fs(xk);
ysk_1   = secant_fs(xk_1);
ysk_new = ysk;
hk      = hk_min;
hk_1    = hk_max;
hk_new  = hk_min;

if (x_history_defined) then
  x_history($+1) = xk
end

if (Log) then
  printf('ls_secant\n');
  yk = secant_f(xk);
  printf('%3s %10s %10s %10s %10s\n', 'It', '|d|', 'x', 'f(x)', 'fs(x)');
  printf('%3d %10.2e %10.2e %10.2e %10.2e\n', 0, 0.0, xk, yk, d'*ysk);
end

for It=1:ItMX
	if (norm(hk-hk_1)<TOL) then
	  if (Log) then
	    printf('ls_secant: convergence level reached\n');
	  end
    hs = hk_new;
    break;  
	end

	hk_new  = hk - (hk - hk_1) / (d'*ysk - d'*ysk_1) * d'*ysk; // Provoque des division par 0 !!
	
	xk_new  = x0 + hk_new*d;
  ysk_new = secant_fs(xk_new);

	xk_1 = xk;     ysk_1 = ysk;
	xk   = xk_new; ysk   = ysk_new;

	hk_1 = hk;
	hk   = hk_new;
	  
  hs = hk_new;
  
  if (x_history_defined) then
    x_history($+1) = xk_new;
  end

	if (Log) then 
	  yk = secant_f(xk);
		printf('%3d %10.2e %10.2e %10.2e %10.2e\n', It, norm(hk_new*d), xk, yk, d'*ysk);
	end
end

// We transform hs as a list for coherency with other methods
aux     = hs;
hs      = list();
hs($+1) = aux;

if (Log) then 
  printf('ls_secant: no convergence in %d iterations\n', ItMX);
end
endfunction

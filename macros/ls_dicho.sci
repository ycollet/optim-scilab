function [hs,x_history] = ls_dicho(dicho_f,dicho_fs,dicho_fss,x0,d,h,Log,TOL,ItMX)
// method of dichotomy for f' 
//
// [hs,x_history]=ls_dicho(f,fs,fss,x0,h,Log,TOL,ItMX)
//
// dicho_f   - handle to function
// dicho_fs  - handle to function evaluating f'
// dicho_fss - handle to function evaluating f'' (ignored)
// x0   - initial guess
// h    - initial interval length
// Log  - print log information, 1=on (default), 0=off
// TOL  - convergence tolerance: (x3-x1)<TOL, default: 1.0e-8
// ItMX - maximal number of iterations, default: 25
// x_history - list of computed points in the parameter space (optional parameter)

StepTOL = 1e-4; // For convergence: we stop the search when the size of the interval is below StepTOL

[nargout, nargin] = argn();

if (nargin < 6) then Log  = 1;                 end
if (nargin < 7) then TOL  = 1.0e-8;            end
if (nargin < 8) then ItMX = 25;                end

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end

xk  = x0;
y0  = dicho_f(xk);
ys0 = dicho_fs(xk);
h1s = h;

if (x_history_defined) then
  x_history($+1) = xk;
end

// First part: localization of the optimum

if (x_history_defined) then
  [h_list_aux, x_history_aux] = ls_2pts_bracket(dicho_f,dicho_fs,dicho_fss,x0,d,h,Log,TOL,ItMX);
  hk_min = h_list_aux(1);
  hk_max = h_list_aux(2);
  x_history = lstcat(x_history, x_history_aux);
else
  [h_list_aux] = ls_2pts_bracket(dicho_f,dicho_fs,dicho_fss,x0,d,h,Log,TOL,ItMX);
  hk_min = h_list_aux(1);
  hk_max = h_list_aux(2);
end

xk_min  = x0 + hk_min*d;
xk_max  = x0 + hk_max*d;
ysk_min = d'*dicho_fs(xk_min);
ysk_max = d'*dicho_fs(xk_max);

if (Log) then
  printf('ls_dicho\n');  
  printf('%3s %10s %15s %15s\n', 'It', 'h', 'd''*fs(x0+h_min*d)', 'd''*fs(x0+h_max*d)');
  printf('%3d %10.2e %15.2e %15.2e\n', 0, 0, ysk_min, ysk_max);
end

if ((hk_min==0)&(hk_max==0)) then
  printf('ls_dicho: warning - initial direction descent not a direction descent\n');
  hs = 0;
  return;
end

// Second part: dichotomy

for It=1:ItMX
  hk_mid  = (hk_max + hk_min) / 2.0;
  xk_mid  = x0 + d*hk_mid;
  ysk_mid = d'*dicho_fs(xk_mid);
  
  if (x_history_defined) then
    x_history($+1) = xk_mid;
  end

  // We perform a dichotomy around the interval we have located
  // We perform a dichotomy of the function f(x)
  
  if (ysk_mid<0) then
    ysk_min = ysk_mid;
    hk_min  = hk_mid;
    xk_min  = xk_mid;
  else
    ysk_max = ysk_mid;
    hk_max  = hk_mid;
    xk_max  = xk_mid;
  end
  
	if (Log) then 
		printf('%3d %10.2e %15.2e %15.2e\n', It, hk_mid, ysk_min, ysk_max);
	end
  	
  hs = hk_mid;
  
  // Test of the validity of the step and test f convergence of the line search
  if (norm(ysk_max-ysk_min)<StepTOL)|(abs(norm(ysk_mid))<TOL) then
    if (Log) then
      printf('ls_dicho: second part terminated after %d iterations\n', It);
    end
    return;
  end
end

// We transform hs as a list (for coherency with other methods)
aux     = hs;
hs      = list();
hs($+1) = aux;

if (Log) then
  printf('ls_dicho: no convergence in %d iterations\n', ItMX);
end
endfunction

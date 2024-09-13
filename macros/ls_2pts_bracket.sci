function [h_list, x_history]=ls_2pts_bracket(ls_2pts_f,ls_2pts_fs,ls_2pts_fss,x0,d,h,Log,TOL,ItMX)
// method of 2 points bracketing for f
//
// [h_list, [x_history]]=ls_2pts_bracket(f,fs,fss,x0,h,Log,TOL,ItMX)
//
// ls_2pts_f   - handle to function
// ls_2pts_fs  - handle to function evaluating f'
// ls_2pts_fss - handle to function evaluating f'' (not used)
// x0   - initial guess
// h    - initial interval length
// Log  - print log information, 1=on (default), 0=off
// TOL  - convergence tolerance: (x3-x1)<TOL, default: 1.0e-8
// ItMX - maximal number of iterations, default: 25
// history - list of computed points in the parameter space
// Output: h_list - a list of step size: h_list(1) = h_min, h_list(2) = h_max
//         x_history - a history of points searched by the method

[nargout, nargin] = argn();

if (nargin < 6) then Log  = 1;                 end
if (nargin < 7) then TOL  = 1.0e-8;            end
if (nargin < 8) then ItMX = 25;                end

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end

// Localization of the optimum

incr = 2.0;

y0  = ls_2pts_f(x0);
ys0 = d'*ls_2pts_fs(x0);

hk_max = h;
hk_min = 0.0;

if (ys0>0) then
  printf('ls_2pts_bracket: warning - initial direction descent not a descent direction\n');
  h_min = 0;
  h_max = 0;
  h_list = list();
  h_list(1) = h_min;
  h_list(2) = h_max;
  return;
end


xk_max = x0 + d*hk_max;
xk_min = x0 + d*hk_min;

if (x_history_defined) then
  x_history = list();
  x_history($+1) = xk_max;
end

ysk_max = d'*ls_2pts_fs(xk_max);
ysk_min = d'*ls_2pts_fs(xk_min);

// Looking for a good starting h_max: a h_max such that fs(x0+h_max*d)<0
if (Log) then 
  printf('ls_2pts_bracket: searching initial h_max\n');
  printf('%3s %10s %10s %10s\n', 'It', 'h', 'f(x0+h*d)', 'd''*fs(x0+h*d)');
  printf('%3d %10.2e %10.2e %10.2e\n', 0, hk_max, ls_2pts_f(xk_max), ysk_max);
end

// Locating the max bound

It_start = 0;

while((ysk_max<0)&(It_start<ItMX))
  It_start = It_start + 1;
  if (x_history_defined) then
    x_history = list();
    x_history($+1) = xk_max;
  end

  if (Log) then 
    printf('%3d %10.2e %10.2e %10.2e\n', It_start, hk_max, ls_2pts_f(xk_max), ysk_max);
  end

  hk_min  = hk_max;
  xk_min  = xk_max;
  ysk_min = ysk_max;

  hk_max  = hk_max * incr;
  xk_max  = x0 + d*hk_max;
  ysk_max = d'*ls_2pts_fs(xk_max);
end

// Sometimes, hk_min can be equal to 0 -> We improve the min boundary for some iterations

if (hk_min==0) then
  It_start = It_start - 2; // We perform at most 2 iterations to improve the min bound
  hk_min = hk_max;
  Begin = %T;
  while(((ysk_min>0) | Begin)&(It_start<ItMX))
    Begin = %F;
    It_start = It_start + 1;
    if (x_history_defined) then
      x_history = list();
      x_history($+1) = xk_min;
    end

    if (Log) then
      printf('%3d %10.2e %10.2e %10.2e\n', It_start, hk_min, ls_2pts_f(xk_min), ysk_min);
    end

    hk_min  = hk_min / incr;
    xk_min  = x0 + d*hk_min;
    ysk_min = d'*ls_2pts_fs(xk_min);
  end
end

if (Log) then 
  printf('ls_2pts_bracket: bracketing\n');
  printf('%3s %10s %10s %10s %10s\n', 'It', 'h_max','h_min','d''*fs(h_max)','d''*fs(h_min)');
  printf('%3d %10.2e %10.2e %10.2e %10.2e\n', It_start, hk_max, hk_min, ysk_max, ysk_min);
end

for It=It_start:ItMX
  // We are looking for a change of sign of the function fs(x0+d*h_max)*fs(x0+d*h_min))
  if (norm(xk_max-xk_min)>TOL) then
    hk_test = (hk_max + hk_min) / 2;
    xk_test = x0 + d*hk_test;
    ysk_test = d'*ls_2pts_fs(xk_test);
    
    if (ysk_test>0) then
      xk_max  = xk_test;
      ysk_max = ysk_test;
      hk_max  = hk_test;
    else
      xk_min  = xk_test;
      ysk_min = ysk_test;
      hk_min  = hk_test;
    end
      
    if (x_history_defined) then
      x_history($+1) = xk_max;
    end
    
    if (Log) then 
      printf('%3d %10.2e %10.2e %10.2e %10.2e\n', It, hk_max, hk_min, ysk_max, ysk_min);
    end
  else
    ItEnd = It;
    break;
  end
  
  ItEnd = It;    
end

h_list = list();
h_list(1) = hk_min;
h_list(2) = hk_max;

if (Log) then
  printf('ls_2pts_bracket: terminated after %d iterations\n', ItEnd);
end
endfunction

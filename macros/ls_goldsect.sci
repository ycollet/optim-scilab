function [hs,x_history] = ls_goldsect(goldsect_f,goldsect_fs,goldsect_fss,x0,d,h,Log,TOL,ItMX)
// method of golden setion dichotomy for f' 
//
// [hs, [x_history]]=ls_goldsect(goldsect_f,goldsect_fs,goldsect_fss,x0,h,Log,TOL,ItMX,cond_func)
//
// goldsect_f   - handle to function
// goldsect_fs  - handle to function evaluating f'
// goldsect_fss - handle to function evaluating f'' (ignored)
// x0   - initial guess
// h    - initial interval length
// Log  - print log information, 1=on (default), 0=off
// TOL  - convergence tolerance: (x3-x1)<TOL, default: 1.0e-8
// ItMX - maximal number of iterations, default: 25
// x_history - list of computed points in the parameter space

F2 = 2/(1+sqrt(5));
F1 = 1 - F2;
StepTOL = 1e-4;

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
y0  = goldsect_f(xk);
ys0 = d'*goldsect_fs(xk);

if (Log) then
  printf('ls_goldsect\n');
  printf('%3s %10s %10s %10s\n', 'It', 'h', 'f_min(h)', 'f_max(h)');
  printf('%3d %10.2e %10.2e %10.2e\n', 0, 0.0, y0, y0);
end

if (x_history_defined) then
  x_history($+1) = x0;
end

// First part: localization of the optimum

if (x_history_defined) then
  [h_list_aux, x_history_aux] = ls_3pts_bracket(goldsect_f,goldsect_fs,goldsect_fss,x0,d,h,Log,TOL,ItMX);
  hk_min = h_list_aux(1);
  hk_mid = h_list_aux(2);
  hk_max = h_list_aux(3);
  x_history = lstcat(x_history, x_history_aux);
else
  [h_list_aux] = ls_3pts_bracket(goldsect_f,goldsect_fs,goldsect_fss,x0,d,h,Log,TOL,ItMX);
  hk_min = h_list_aux(1);
  hk_mid = h_list_aux(2);
  hk_max = h_list_aux(3);
end

if ((hk_min==0)&(hk_max==0)&(hk_mid==0)) then
  printf('ls_goldsect: warning - initial direction descent not a direction descent\n');
  hs = 0;
  return;
end

xk_min = x0 + hk_min*d;
xk_max = x0 + hk_max*d;
yk_min = goldsect_f(xk_min);
yk_max = goldsect_f(xk_max);

// Second part: golden section dichotomy

if (Log) then
  printf('ls_goldsect\n');
  printf('%3s %10s %10s %10s\n', 'It', 'h', 'f_min(h)', 'f_max(h)');
end

for It=1:ItMX
  hk1_test = hk_min + F1*(hk_max - hk_min);
  hk2_test = hk_min + F2*(hk_max - hk_min);
  
  xk1_test = x0+hk1_test*d;
  xk2_test = x0+hk2_test*d;

  yk1_test = goldsect_f(xk1_test);
  yk2_test = goldsect_f(xk2_test);
  
  if (yk1_test<yk2_test) then
    xk_max = xk2_test;
    yk_max = yk2_test;
    hk_max = hk2_test;
  elseif (yk1_test>yk2_test) then
    xk_min = xk1_test;
    yk_min = yk1_test;
    hk_min = hk1_test;
  else
    xk_max = xk2_test;
    yk_max = yk2_test;
    hk_max = hk2_test;

    xk_min = xk1_test;
    yk_min = yk1_test;
    hk_min = hk1_test;
  end
  
  if (Log) then 
    printf('%3d %10.2e %10.2e %10.2e %10.2e %10.2e\n', It, hk_min, yk_min, yk_max, yk1_test, yk2_test);
  end

  if (x_history_defined) then
    x_history($+1) = xk1_test;
    x_history($+1) = xk2_test;
  end
    
  // Computation of the step size
  hs = hk_min;
  
  // Test of the validity of the step
  if (abs(norm(goldsect_fs(xk_min)))<TOL)|(abs(norm(xk_max-xk_min))<StepTOL) then
    if (Log) then
      printf('ls_goldsect: second part terminated after %d iterations\n', It);
    end
    return;
  end
end

// We transform hs as a list for coherency with other methods
aux     = hs;
hs      = list();
hs($+1) = aux;

if (Log) then
  printf('ls_goldsect: no convergence in %d iterations\n', ItMX);
end
endfunction

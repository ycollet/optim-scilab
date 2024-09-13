function [hs, x_history]=ls_backtrack(backtrack_f,backtrack_fs,backtrack_fss,x0,d,h,Log,TOL,ItMX);
// method of backtracking for f' 
//
// [hs,[x_history]]=ls_backtrack(backtrack_f,backtrack_fs,backtrack_fss,x0,h,Log,TOL,ItMX)
//
// backtrack_f   - handle to function
// backtrack_fs  - handle to function evaluating f'
// backtrack_fss - handle to function evaluating f'' (ignored)
// x0   - initial guess
// h    - initial interval length
// Log  - print log information, 1=on (default), 0=off
// TOL  - convergence tolerance: (x3-x1)<TOL, default: 1.0e-8
// ItMX - maximal number of iterations, default: 25
// x_history - list of computed points in the parameter space

[nargout, nargin] = argn();

if ~isdef('Log','local')       then Log  = 1;                 end
if ~isdef('TOL','local')       then TOL  = 1.0e-8;            end
if ~isdef('ItMX','local')      then ItMX = 25;                end

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end

cg_alpha = 0.0001;
cg_beta  = 0.8;
c1       = 1E-4; // c1 in [0,1]

y0  = backtrack_f(x0);
ys0 = backtrack_fs(x0);

if (Log) then
  printf('ls_backtrack\n');
  printf('%3s %10s %10s %10s %10s\n', 'It', '|d|', 'x', 'f(x)', 'fs(x)');
  printf('%3d %10.2e %10.2e %10.2e %10.2e\n', 0, 0.0, x0, y0, d'*ys0);
end

// backtracking algorithm

// Initialization of t
//t = - ys0' * d / norm(d)^2;
t = 2;
if (Log) then
  printf('ls_backtrack: initial value of t = %f\n', t);
end

Index = 1;

if (x_history_defined) then
  x_history($+1) = x0+t*d;
end

// While the step doesn't verify the armijo condition, we reduce the step
yk  = backtrack_f(x0+t*d);
while ((Index<ItMX) & (yk <= y0 + c1*t*ys0'*d))
  t = t * cg_beta;
  Index = Index + 1;
  
  if (Log) then
    printf('%3d %10.2e %10.2e %10.2e %10.2e\n', Index, t, x0 - t/cg_beta * d, yk, d'*yk);
  end
  
  yk  = backtrack_f(x0+t*d);

  if (x_history_defined) then
    x_history($+1) = x0+t*d;
  end
end

hs      = list();
hs($+1) = t;
endfunction

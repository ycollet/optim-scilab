function [x_opt, x_history] = optim_dfp(x0, dfp_f, dfp_df, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart)
//
//////////////////
// Informations //
//////////////////
//
// DFP optimization method
// 
// x0       : initial starting point (must be a column vector)
// dfp_f    : objective function
// dfp_df   : gradient of the objective function (must return a column vector)
// h        : size of the initial step (optional parameter)
// Log      : verbose mode or not (%T / %F) (optional parameter)
// ItMX     : number of steepest descent step (optional parameter)
// ls_ItMX  : number of iterations of the line search (optional parameter)
// lsm      : type of line search method (optional parameter)
// TOL      : accuracy for convergence test - derivatives (optional parameter)
// StepTOL  : accuracy for convergence test - size of the step (optional parameter)
// XTOL     : accuracy for convergence test - improvement on x between two iterations (optional parameter)
// SubTOL   : accuracy for convergence test (line search) (optional parameter)
// Restart  : after how many iterations do we restart the search direction

///////////////////////////
// Testing the arguments //
///////////////////////////

[nargout, nargin] = argn();

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history         = list();
  x_history($+1)    = list();
  x_history($)($+1) = x0;
end

if ~isdef('x0','local') then
  error('optim_dfp: x0 is mandatory');
end
if (~isdef('h','local')) then
  h = 0.00125; // initial step length in line search
end
if (~isdef('Log','local')) then
  Log = %F; // print info from line step, %T = on / %F = off
end
if (~isdef('ItMX','local')) then
  ItMX = 50; // maximum number of descent steps
end
if (~isdef('ls_ItMX','local')) then
  ls_ItMX = 100; // maximum number of line search iteration  
end
if (~isdef('lsm','local')) then
  lsm = ls_goldsect; // type of line search method used
end
if (~isdef('TOL','local')) then
  TOL = 1.0e-4; // accuracy for convergence test - derivatives 
end
if (~isdef('StepTOL','local')) then
  StepTOL = 1.0e-4; // accuracy for convergence test - size of the step   
end
if (~isdef('XTOL','local')) then
  XTOL = 1.0e-4; // accuracy for convergence test - improvement on x between two iterations
end
if (~isdef('SubTOL','local')) then
  SubTOL = 1.0e-2; // accuracy for convergence test (line search)
end

if ~isdef('dfp_f','local') then
  error('optim_dfp: dfp_f is mandatory');
else
  if typeof(dfp_f)=='list' then
    deff('y=_dfp_f(x)','y=dfp_f(1)(x, dfp_f(2:$))');
  else
    deff('y=_dfp_f(x)','y=dfp_f(x)');
  end
end
if ~isdef('dfp_df','local') then
  error('optim_dfp: dfp_df is mandatory');
else
  if typeof(dfp_df)=='list' then
    deff('y=_dfp_df(x)','y=dfp_df(1)(x, dfp_df(2:$))');
  else
    deff('y=_dfp_df(x)','y=dfp_df(x)');
  end
end

//////////////////////////////////////////
// Davidon Fletcher Powell quasi-Newton //
//////////////////////////////////////////

n = size(x0,1);
step = [];
yk_1 = [];
xk_1 = x0;
yk_2 = [];
xk_2 = x0;
dk   = [];
dyk  = [];
H    = [];

for It=1:ItMX
  if (Log) then
    printf('DFP Iteration %d - ', It);
  end
  if ((modulo(It, Restart)==0)|(It==1)) then
    if (Log) then
      printf('DFP restart :\n');
    end
    H    = eye(n,n);
    yk_1 = _dfp_df(xk_1);
    d    = -H*yk_1;
  else
    if (Log) then
      printf('DFP iteration :\n');
    end
    yk_2 = yk_1;
    yk_1 = _dfp_df(xk_1);
    dk   = xk_1 - xk_2;
    dyk  = yk_1 - yk_2;

    H = H + (dk*dk')/(dk'*dyk) - H*dyk*dyk'*H/(dyk'*H*dyk);
    d = -H*yk_1; 
  end

  //
  // Computation of the step size	
  //

  // d = d/norm(d);
  [step, x_history_aux] = lsm(_dfp_f, _dfp_df, H, xk_1, d, h, Log, SubTOL, ls_ItMX);

  if (Log) then
    scf(); plot_d(_dfp_f, xk_1, d, lbounds, ubounds, %T, 0, max([h abs(step(1))]));
    printf('optim_cg: close the graph to continue - Iteration %d / %d\n', It, ItMX);
    xclick();
    xdel();
  end

  if (Log) then
    printf('descent direction:');disp(d');
    printf('step size = %f\n', step(1));
  end

  xk_2 = xk_1;
  xk_1 = xk_1 + step(1)*d;
  
  if (x_history_defined) then
    x_history($+1) = list();
    x_history($)   = lstcat(x_history($), x_history_aux);
  end
  
  ys = _dfp_df(xk_1);
  y  = _dfp_f(xk_1);
	if ((abs(norm(ys))<TOL)) then
	  if Log then
      printf('optim_dfp: convergence reached:\n')
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, norm(ys));
    end
    x_opt = xk_1;
    return;
  elseif ((step(1)<StepTOL)|(norm(xk_1-xk_2)<XTOL)) then
    if Log then
      printf('optim_dfp: optimization stopped:\n')
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, norm(ys));
    end
    x_opt = xk_1;
    return;
  else
    if Log then
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, norm(ys));
    end
  end  
end

x_opt = xk_1; // We return the last point
endfunction

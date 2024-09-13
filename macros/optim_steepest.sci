function [x_opt, x_history] = optim_steepest(x0, steepest_f, steepest_df, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL)
//
//////////////////
// Informations //
//////////////////
//
// Steepest descent optimization method
// 
// x0       : initial starting point (must be a column vector)
// f        : objective function
// df       : gradient of the objective function (must return a column vector)
// h        : size of the initial step (optional parameter)
// Log      : verbose mode or not (%T / %F) (optional parameter)
// ItMX     : number of steepest descent step (optional parameter)
// ls_ItMX  : number of iterations of the line search (optional parameter)
// lsm      : type of line search method (optional parameter)
// TOL      : accuracy for convergence test - derivatives (optional parameter)
// StepTOL  : accuracy for convergence test - size of the step (optional parameter)
// XTOL     : accuracy for convergence test - improvement on x between two iterations (optional parameter)
// SubTOL   : accuracy for convergence test (line search) (optional parameter)
//

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
  error('optim_steepest: x0 is mandatory');
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
if ~isdef('steepest_f','local') then
  error('optim_steepest: steepest_f is mandatory');
else
  if typeof(steepest_f)=='list' then
    deff('y=_steepest_f(x)','y=steepest_f(1)(x, steepest_f(2:$))');
  else
    deff('y=_steepest_f(x)','y=steepest_f(x)');
  end
end
if ~isdef('steepest_df','local') then
  error('optim_steepest: steepest_df is mandatory');
else
  if typeof(steepest_df)=='list' then
    deff('y=_steepest_df(x)','y=steepest_df(1)(x, steepest_df(2:$))');
  else
    deff('y=_steepest_df(x)','y=steepest_df(x)');
  end
end

//////////////////////
// steepest descent //
//////////////////////

xk_1 = x0;
xk_2 = x0;

for It=1:ItMX
  //
  // Computation of the descent direction
  //
  d = -_steepest_df(xk_1);

  //
  // Computation of the step size	
  //
  
  //d = d/norm(d);
  
  [step, x_history_aux] = lsm(_steepest_f, _steepest_df, [], xk_1, d, h, Log, SubTOL, ls_ItMX);

  if (Log) then
    scf(); plot_d(_steepest_f, xk_1, d, lbounds, ubounds, %T, 0, max([h abs(step(1))]));
    printf('optim_steepest: close the graph to continue - Iteration %d / %d\n', It, ItMX);
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
  
  ys = _steepest_df(xk_1);
  y  = _steepest_f(xk_1);
  
	if (norm(ys)<TOL) then
	  if Log then
      printf('optim_steepest: minimum found:\n')
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, norm(ys));
    end
    x_opt = xk_1;
    return;
	elseif ((step(1)<StepTOL)|(norm(xk_1-xk_2)<XTOL)) then
	  if Log then
      printf('optim_steepest: optimization stopped:\n')
      printf('%10s %10s %16s\n', 'f(x)', '||fs(x)||');
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

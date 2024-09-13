function [x_opt, x_history] = optim_cg(x0, cg_f, cg_df, h, Log, ItMX, ls_ItMX, lsm, typecg, TOL, StepTOL, XTOL, SubTOL, Restart, Spectral)

//////////////////
// Informations //
//////////////////

// Steepest descent optimization method
// 
// x0       : initial starting point (must be a column vector)
// cg_f     : objective function
// cg_df    : gradient of the objective function (must return a column vector)
// h        : size of the initial step (optional parameter)
// Log      : verbose mode or not (%T / %F) (optional parameter)
// ItMX     : number of steepest descent step (optional parameter)
// ls_ItMX  : number of iterations of the line search (optional parameter)
// lsm      : type of line search method (optional parameter)
// typecg   : type of conjugate gradient used (optional parameter)
//            possible values are : - 'fr' for Fletcher-Reeves (default value)
//                                  - 'pr' for Polak-Ribiere
//                                  - 'hs' for Hestenes-Siefel
//                                  - 'pw' for Powel
//                                  - 'dy' for Dai-Yuan
// TOL      : accuracy for convergence test - derivatives (optional parameter)
// StepTOL  : accuracy for convergence test - size of the step (optional parameter)
// XTOL     : accuracy for convergence test - improvement on x between two iterations (optional parameter)
// SubTOL   : accuracy for convergence test (line search) (optional parameter)
// Restart  : after how many iterations do we reinitialize the search direction
// Spectral : flag to compute the "spectral correction" ofthe search direction

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

ddf = [];

if ~isdef('x0','local') then
  error('optim_cg: x0 is mandatory');
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
if (~isdef('typecg','local')) then
  typecg = 'fr'; // type of conjugate gradient method used. Default value: fletcher-reeves
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
if (~isdef('Restart','local')) then
  Restart = 7;
end
if (~isdef('Spectral','local')) then
  Spectral = %F;
end

if ~isdef('cg_f','local') then
  error('optim_cg: cg_f is mandatory');
else
  if typeof(cg_f)=='list' then
    deff('y=_cg_f(x)','y=cg_f(1)(x, cg_f(2:$))');
  else
    deff('y=_cg_f(x)','y=cg_f(x)');
  end
end
if ~isdef('cg_df','local') then
  error('optim_cg: cg_df is mandatory');
else
  if typeof(cg_df)=='list' then
    deff('y=_cg_df(x)','y=cg_df(1)(x, cg_df(2:$))');
  else
    deff('y=_cg_df(x)','y=cg_df(x)');
  end
end

nu = 1e-4;

/////////////////////////
// conjugate gradients //
/////////////////////////

g0    = [];
omega = [];
d0    = [];

xk_1 = x0;
xk_2 = x0;

for It=1:ItMX
  if ((modulo(It, Restart)==0)|(It==1)|PowRestart) then
    if (Log) then
      printf('cgdemo: gradient reinitialization \n');
    end
    g1 = _cg_df(xk_1);
    d1 = -g1;
    PowRestart = %F;
  else
    if (Log) then
      printf('cgdemo: conjugate gradient step \n');
    end
    d0 = d1;
    g0 = g1;
    g1 = _cg_df(xk_1);
    select cgm
      case 'fr' then
        // Fletcher-Reeves
        omega = (g1'*g1)/(g0'*g0);
      case 'pr' then
        // Polak-Ribiere
        omega = g1'*(g1 - g0)/(d0'*(g1 - g0));
      case 'hs' then
        // Hestenes-Stiefel
        omega = max([0 g1'*(g1 - g0)/(g1'*g0)]);
      case 'pw' then
        // Powel
        omega = g1'*(g1 - g0)/(g0'*g0);
      case 'dy' then
        // Dai-Yuan
        omega = (g1'*g1)/(g1'*g0);
      else
        error('unknown function %s', cgm);
    end
    
    if (Spectral) then
       sk = xk_1 - xk_2;
       yk = g1 - g0;
      theta = sk'*sk / (sk'*yk);
    else
      theta = 1;
    end
        
    d1 = - theta * g1 + omega * d0;

    // If the search direction is not sufficiently conjugate, we start a new
    // gradient iteration
    PowRestart = ((g1'*g0)/(g0'*g0)>nu);
    // If the step is too small, we start a new gradient iteration
    //PowRestart = PowRestart | ((norm(g1-g0)/norm(g0))<StepIncrease);
  end
	
  //
  // Computation of the step size
  //
  
  //d1 = d1/norm(d1);
  
 	[step, x_history_aux] = lsm(_cg_f, _cg_df, [], xk_1, d1, h, Log, SubTOL, ls_ItMX);

//  if (Log) then
//    scf(); plot_d(_cg_f, xk_1, d1, lbounds, ubounds, %T, 0, max([h abs(step(1))]));
//    printf('optim_cg: close the graph to continue - Iteration %d / %d\n', It, ItMX);
//    xclick();
//    xdel();
//  end

  if (Log) then
    printf('descent direction:');disp(d1');
    printf('step size = %f\n', step(1));
  end

  xk_2 = xk_1;
  if (cgm=='pr') then
    xk_1 = xk_1 + max([0 step(1)]) * d1;
  else
    xk_1 = xk_1 + step(1) * d1;
  end
		
	if (x_history_defined) then
	  x_history($+1) = list();
	  x_history($) = lstcat(x_history($), x_history_aux);
	end
	
	ys = _cg_df(xk_1);
	y  = _cg_f(xk_1);
	
	if (norm(ys)<TOL) then
	 if Log then
      printf('optim_cg: convergence reached:\n')
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, ys);
    end
    x_opt = xk_1;
    return;
	elseif ((step(1)<StepTOL)|norm(xk_1-xk_2)<XTOL) then
	 if Log then
      printf('optim_cg: optimization stopped:\n')
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, ys);
    end
    x_opt = xk_1;
    return;
	else
	  if Log then
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, ys);
    end
  end
end

x_opt = xk_1 // We return the last solution
endfunction

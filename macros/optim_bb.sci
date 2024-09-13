function [x_opt, x_history] = optim_bb(x0, bb_f, bb_df, Log, ItMX, DXMax, TOL, XTOL)
//
//////////////////
// Informations //
//////////////////
//
// Barzilai-Borwein optimization method
// 
// x0       : initial starting point (must be a column vector)
// f        : objective function
// df       : gradient of the objective function (must return a column vector)
// h        : size of the initial step (optional parameter)
// Log      : verbose mode or not (%T / %F) (optional parameter)
// ItMX     : number of step (optional parameter)
// TOL      : accuracy for convergence test - derivatives (optional parameter)
// XTOL     : accuracy for convergence test - improvement on x between two iterations (optional parameter)
//

///////////////////////////
// Testing the arguments //
///////////////////////////

[nargout, nargin] = argn();

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history      = list();
  x_history($+1) = x0;
end

if ~isdef('x0','local') then
  error('optim_bb: x0 is mandatory');
end
if (~isdef('Log','local')) then
  Log = %F; // print info from line step, %T = on / %F = off
end
if (~isdef('ItMX','local')) then
  ItMX = 50; // maximum number of descent steps
end
if (~isdef('TOL','local')) then
  TOL = 1.0e-4; // accuracy for convergence test - derivatives 
end
if (~isdef('XTOL','local')) then
  XTOL = 1.0e-4; // accuracy for convergence test - improvement on x between two iterations
end
if (~isdef('DXMax','local')) then
  DXMax = 1;
end
if ~isdef('bb_f','local') then
  error('optim_bb: bb_f is mandatory');
else
  if typeof(bb_f)=='list' then
    deff('y=_bb_f(x)','y=bb_f(1)(x, bb_f(2:$))');
  else
    deff('y=_bb_f(x)','y=bb_f(x)');
  end
end
if ~isdef('bb_df','local') then
  error('optim_bb: bb_df is mandatory');
else
  if typeof(bb_df)=='list' then
    deff('y=_bb_df(x)','y=bb_df(1)(x, bb_df(2:$))');
  else
    deff('y=_bb_df(x)','y=bb_df(x)');
  end
end

//////////////////////
// Barzilai-Borwein //
//////////////////////

xk_1 = x0;
xk_2 = x0;

alpha = 1;

for It=1:ItMX
  //
  // Computation of the descent direction
  //
  
  d = -_bb_df(xk_1)/alpha;
      
  xk_2 = xk_1;
  
  if (norm(d)<=DXMax) then
  	 xk_1 = xk_1 + d;
  else 
    xk_1 = xk_1 + d/(norm(d)) * DXMax;
  end
	 
  yk    = (_bb_df(xk_1) - d*alpha);
  alpha = abs((yk' * yk) / (yk' * d));
 
  if (x_history_defined) then
    x_history($+1) = xk_1;
  end
  
  if Log then
    ys = _bb_df(xk_1);
    y  = _bb_f(xk_1);
  end

  if Log & It==1 then
    printf('%10s %10s %10s\n', 'f(x)', '||fs(x)||','alpha');
  else
    printf('%10.2e %10.2e %10.2e\n', y, norm(ys),alpha);
  end
  
	if (norm(ys)<TOL) then
	  if Log then
      printf('optim_bb: derivative convergence criterion reached\n')
    end
    x_opt = xk_1;
    return;
	end
	if (norm(xk_2 - xk_1)<XTOL) then
	  if Log then
	    printf('optim_bb: displacement convergence criterion reached\n');
	  end
	  x_opt = xk_1;
	  return;
	end
end

x_opt = xk_1; // We return the last point
endfunction

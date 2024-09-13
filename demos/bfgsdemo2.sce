function [x_opt, x_history, H_history] = optim_bfgs(x0, bfgs_f, bfgs_df, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart)

[nargout, nargin] = argn();

H_history_defined = (nargout==3);
x_history_defined = (nargout>=2);

if (x_history_defined) then
  x_history         = list();
  x_history($+1)    = list();
  x_history($)($+1) = x0;
end

if H_history_defined then
  H_history = list();
end

if ~isdef('bfgs_f','local') then
  error('optim_bfgs: bfgs_f is mandatory');
else
  if typeof(bfgs_f)=='list' then
    deff('y=_bfgs_f(x)','y=bfgs_f(1)(x, bfgs_f(2:$))');
  else
    deff('y=_bfgs_f(x)','y=bfgs_f(x)');
  end
end
if ~isdef('bfgs_df','local') then
  error('optim_bfgs: bfgs_df is mandatory');
else
  if typeof(bfgs_df)=='list' then
    deff('y=_bfgs_df(x)','y=bfgs_df(1)(x, bfgs_df(2:$))');
  else
    deff('y=_bfgs_df(x)','y=bfgs_df(x)');
  end
end

///////////////////////////////////////////////////
// Broyden Fletcher Goldfarb Shanno quasi-Newton //
///////////////////////////////////////////////////

xk_1 = x0;
xk_2 = x0;

n  = size(x0,1);

// Most robust initialisation for several cases
// Doesn't seem to be the best one for sphere function.
// Add some random to get the best estimate
In = eye(n,n);
//In = diag(0.01*rand(1,n)) + eye(n,n);
//In = diag(rand(1,n));
//In = diag((1-0.3)*rand(1,n)+0.3);

for It=1:ItMX
  //
  // Computation of the descent direction
  //
  if (Log) then
    printf('BFGS Iteration %d - ', It);
  end
  if ((modulo(It, Restart)==0)|(It==1)) then
    H    = In;
    yk_1 = _bfgs_df(xk_1);
    d    = -H*yk_1;
  else
    if (Log) then
      printf('BFGS iteration :\n');
    end
    yk_2 = yk_1;
    yk_1 = _bfgs_df(xk_1);
    dk   = xk_1 - xk_2;
    dyk  = yk_1 - yk_2;

//    H = H + (dyk*dyk')/(dyk'*dk) - (H*dk*dk'*H')/(dk'*H*dk);
//    d = -H^(-1)*yk_1; // No the best way. We should apply the BFGS relationship which gives
//                      // the inverse of the H matrix
    // Sherman–Morrison formula
    H = H + (dk*dk')*(dk'*dyk+dyk'*H*dyk)/(dk'*dyk)^2-(H*dyk*dk'+dk*dyk'*H)/(dk'*dyk);
    d = -H*yk_1; // No the best way. We should apply the BFGS relationship which gives
                      // the inverse of the H matrix

  end
  
  if H_history_defined then
    H_history($+1) = H;
  end
    
  //
  // Computation of the step size	
  //
  
  // d = d / norm(d);
  
  [step, x_history_aux] = lsm(_bfgs_f, _bfgs_df, H, xk_1, d, h, Log, SubTOL, ls_ItMX);

  // Computation of the step
  xk_2 = xk_1;
  xk_1 = xk_1 + step(1)*d;
  
  if (x_history_defined) then
    x_history($+1) = list();
    x_history($)   = lstcat(x_history($), x_history_aux);
  end
  
  ys = _bfgs_df(xk_1);
  y  = _bfgs_f(xk_1);
	if ((abs(norm(ys))<TOL)) then
    x_opt = xk_1;
    return;
	elseif ((step(1)<StepTOL)|(abs(norm(xk_1-xk_2))<XTOL)) then
    x_opt = xk_1
    return;
  end  
end

x_opt = xk_1; // We return the last point
endfunction

//////////////////
// Main program //
//////////////////


deff('y=f(x)','y=gen_rosen(x,10)');
x0 = ones(10,1); x0(1) = -1;
//deff('y=f(x)','y=sphere(x,10)');
//x0 = 4*rand(10,1);
deff('y=df(x)','y=derivative(f,x)''');

ItMX     = 100;         // maximum number of descent steps
ls_ItMX  = 100;        // maximum number of line search iteration
Restart  = ItMX;       // restart cg every Restart iterations
TOL      = 1e-6;     // accuracy for convergence test (minimum)
StepTOL  = 1e-12;    // accuracy for convergence test - size of the step 
XTOL     = 1e-12;    // accuracy for convergence test - improvement on x between two iterations
SubTOL   = 1e-6;     // accuracy for convergence test (line search)
Log      = %F;         // print info from line step, 1 = on / 0 = off
h        = 4.0; //0.01;       // initial step length in line search
// lsm  = ls_secant;
lsm  = ls_dicho;
// lsm  = ls_goldsect;
// lsm  = ls_backtrack;
// lsm  = ls_polynom;

////////////////////////////////////////////////////////////////////////////////////////////////////

////////////// Broyden Fletcher Goldfarb Shanno quasi-Newton

[x_opt, x_history, H_history] = optim_bfgs(x0, f, df, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL);

t = []; Y = [];
t = 1:length(x_history)-1;
for i=1:length(x_history)-1
  [tmp,H_diff] = derivative(f,x_history(i)(1),H_form='hypermat');
  //Y(i) = norm(H_diff - H_history(i));
  Y(i) = norm(H_diff^(-1) - H_history(i));
end

t_diff = []; Y_diff = [];
t_diff = 1:length(x_history)-2;
for i=2:length(x_history)-1
  [tmp,H_diff] = derivative(f,x_history(i)(1),H_form='hypermat');
  Y_diff(i-1) = norm(H_history(i-1) - H_history(i));
end

subplot(2,1,1);
plot(t,Y,'k-');
xtitle('Difference between real Hessian and estimated one','Iteration','norm');

subplot(2,1,2);
plot(t_diff,Y_diff,'k-');
xtitle('Progress in Hessian estimation','Iteration','norm');

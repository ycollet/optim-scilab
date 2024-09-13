function [x_opt, x_history] = optim_lbfgs(x0, bfgs_f, bfgs_df, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart, m)
//
//////////////////
// Informations //
//////////////////
//
// BFGS optimization method
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
// Restart  : after how many iteration do we restart the search direction
// m        : number of history steps

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
  error('optim_bfgs: x0 is mandatory');
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
if (~isdef('Restart','local')) then
  Restart = 7;
end
if (~isdef('m','local')) then
  m = 4;
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

/////////////////////////
// Limited memory BFGS //
/////////////////////////

xk   = x0 ;
xk_1 = x0 ;
k    = 1;

n   = size(x0,1);
In  = eye(n,n);
Y   = [];
S   = [];
rho = [];

for It=1:ItMX
  //
  // Computation of the descent direction
  //
  // LBFGS
  
  if (Log) then
    printf('BFGS Iteration %d - ', It);
  end

  if ((modulo(It, Restart)==0)|(It==1)) then
    if (Log) then
      printf('BFGS restart :\n');
    end
    H  = In;
    gk = _bfgs_df(xk); 
    d  = -gk;
  else
    if (Log) then
      printf('BFGS iteration :\n');
    end
    H      = In;
    gk_1   = _bfgs_df(xk_1); 
    d      = -H*gk_1;
    gk     =  _bfgs_df(xk);
    sk     = xk - xk_1;
    yk     = gk - gk_1;
    vm_rho = 1/(yk'  *sk);
    Vk     = In - vm_rho*yk*sk';
    mc     = min(It-1,m);

    if (It-1<m) then
      if isempty(Y) then
        Y(:,1) = yk;
      else
        Y(:,$+1) = yk;
      end
    else
      Y(:,1:$-1) = Y(:,2:$);
      Y(:,$) = yk;
    end
    
    if (It-1<m) then
      if isempty(S) then
        S(:,1) = sk;
      else
        S(:,$+1) = sk;
      end
    else
      S(:,1:$-1) = S(:,2:$);
      S(:,$) = sk;
    end

    if (It-1<m) then
      if isempty(rho) then
        rho(1) = vm_rho;
      else
        rho($+1) = vm_rho;
      end
    else
      rho(1:$-1) = rho(2:$);
      rho($) = vm_rho;
    end

    Vect_tmp = In;
    
    n_vec = size(S,2);
    
    for i=1:n_vec
      Vect_tmp = Vect_tmp * (In - rho(i) * Y(:,i) * S(:,i)');
    end
    H_tmp = Vect_tmp' * H * Vect_tmp;
       
    for i=1:n_vec
      Vect_tmp = In;
      for j=1:n_vec-i+1
        Vect_tmp = Vect_tmp * (In - rho(n_vec-j+1) * Y(:,n_vec-j+1) * S(:,n_vec-j+1)');
      end
      H_tmp = H_tmp + rho(n_vec+1-i)*Vect_tmp'*S(:,n_vec-i+1)*S(:,n_vec-i+1)'*Vect_tmp
    end
    
    H_tmp = H_tmp + rho(n_vec)*S(:,n_vec)*S(:,n_vec)'

    gk = _bfgs_df(xk); 
    H  = H_tmp;
    d  = -H'*gk;
  end
  
  //
  // Computation of the step size	
  //
  
  [step, x_history_aux] = lsm(_bfgs_f, _bfgs_df, H, xk, d, h, Log, SubTOL, ls_ItMX);

  if (Log) then
    printf('descent direction:');disp(d');
    printf('step size = %f\n', step(1));
  end

  // Computation of the step
  xk_1 = xk;
  xk   = xk + step(1)*d;
  
  if (x_history_defined) then
    x_history($+1) = list();
    x_history($)   = lstcat(x_history($), x_history_aux);
  end
  
  gk = _bfgs_df(xk);
  y  = _bfgs_f(xk);
	if ((abs(norm(gk))<TOL)) then
	  if Log then
      printf('\nminimum found:\n')
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, norm(gk));
    end
    x_opt = xk;
    return;
	elseif ((step(1)<StepTOL)|(abs(norm(xk-xk_1))<XTOL)) then
	  if Log then
      printf('\noptimization stopped:\n')
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, norm(gk));
    end
    x_opt = xk
    return;
  else
    if Log then
      printf('%10s %10s\n', 'f(x)', '||fs(x)||');
      printf('%10.2e %10.2e\n', y, norm(gk));
    end
  end  
end

x_opt = xk; // We return the last point
endfunction

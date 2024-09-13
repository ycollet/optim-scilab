function [x_opt, x_history] = optim_feasdir(fdir_f, fdir_df, fdir_g, fdir_dg, x0, ...
                                            ItMX, MaxEvalFunc, upper, lower, Theta_0, epsilon, ...
                                            h, MaxMinStepCount, SubTOL, ls_ItMX, Log)

if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('MaxMinStepCount','local') then
  MaxMinStepCount = 10;
end
if (~isdef('h','local')) then
  h = 0.00125; // initial step length in line search
end
if (~isdef('ls_ItMX','local')) then
  ls_ItMX = 100; // maximum number of line search iteration  
end
if (~isdef('lsm','local')) then
  lsm = ls_goldsect; // type of line search method used
end
if (~isdef('SubTOL','local')) then
  SubTOL = 1.0e-2; // accuracy for convergence test (line search)
end
if (~isdef('Theta_0','local')) then
  Theta_0 = 1.0;
end
if (~isdef('epsilon','local')) then
  epsilon = -0.01;
end
if ~isdef('fdir_f','local') then
  error('optim_feasdir: fdir_f is mandatory');
else
  if typeof(fdir_f)=='list' then
    deff('y=_fdir_f(x)','y=fdir_f(1)(x, fdir_f(2:$))');
  else
    deff('y=_fdir_f(x)','y=fdir_f(x)');
  end
end
if ~isdef('fdir_df','local') then
  error('optim_feasdir: fdir_df is mandatory');
else
  if typeof(fdir_df)=='list' then
    deff('y=_fdir_df(x)','y=fdir_df(1)(x, fdir_df(2:$))');
  else
    deff('y=_fdir_df(x)','y=fdir_df(x)');
  end
end

if ~isdef('fdir_g','local') then
  error('optim_feasdir: fdir_g is mandatory');
else
  if typeof(fdir_g)=='list' then
    deff('y=_fdir_g(x)','y=fdir_g(1)(x, fdir_g(2:$))');
  else
    deff('y=_fdir_g(x)','y=fdir_g(x)');
  end
end

if ~isdef('fdir_dg','local') then
  error('optim_feasdir: fdir_dg is mandatory');
else
  if typeof(fdir_dg)=='list' then
    deff('y=_fdir_dg(x)','y=fdir_dg(1)(x, fdir_vdg(2:$))');
  else
    deff('y=_fdir_dg(x)','y=fdir_dg(x)');
  end
end

g_value  = _fdir_g(x0);
dg_value = _fdir_dg(x0);

NbIneqConstr = max(size(g_value));

[nargout, nargin] = argn();

x_history_defined = (nargout>=2);
if (x_history_defined) then
  x_history = list();
  x_history($+1) = list();
  x_history($)($+1) = x0;
end

// dg must return a matrix NbVar x NbConstr
// x0 must be a column vector of size NbVar
// g must be a column vector of size NbConstr
// df must be a column vector of size NbVar
// delta_ml can either be a scalar or a column vector of size NbVar
// lower and upper must be column vectors of size NbVar

IneqConstrMatr = [];
IneqBoundMatr  = [];

f_value  = _fdir_f(x0);
df_value = _fdir_df(x0);
NbVar    = max(size(x0)); 

x        = x0;
x_k_m_1 = x0;

Iteration  = 0;
EvalFunc   = 0;
XTol_Count = 0;
XTol       = %inf;

Debug = %F;

// Verification of the size of all the parameters
if (size(x0,1)~=NbVar) then
  error('x0 must be a column vector of size nbvar');
end

if size(g_value,1)~=NbIneqConstr then
  error('g must return a column vector of size nbconstr');
end
if ((size(dg_value,1)~=NbVar)&(size(dg_value,2)~=NbIneqConstr)) then
  error('dg must return a matrix of size nbvar x nbconstr');
end

if (size(df_value,1)~=NbVar) then
  error('df must return a column vector of size nbvar');
end
if (size(lower,1)~=NbVar) then
  error('lower must be a column vector of size NbVar');
end
if (size(upper,1)~=NbVar) then
  error('upper must be a column vector of size NbVar');
end

if or(g_value>0) then
  if (Log) then
    warning('optim_feasdir: Initial solution not feasible. Find a new feasible starting point.');
  end
  xopt = x;
  return
end

Iteration = 0;

while((Iteration<ItMX)&(EvalFunc<MaxEvalFunc))
  Iteration = Iteration + 1;
  EvalFunc  = EvalFunc + 1;
  x_k_m_1   = x;
  
  if (x_history_defined) then
    x_history($+1) = list();
    x_history($)($+1) = x;
  end

  // Building local linear models
  // For the objective function
  if (Iteration~=1) then
    f_value  = _fdir_f(x);
    df_value = _fdir_df(x);
    g_value  = _fdir_g(x);
    dg_value = _fdir_dg(x);
  end
  // Objective function matrix
  ObjFuncMatr = [zeros(1,length(x)) -1]'; // We maximize the objective function
  
  // Building the matrix constraint
  // The Inequality constraints
  Index_ActivConstr = find(g_value>0);
  Theta = (1 - g_value/epsilon).^2 * Theta_0;
  //Theta = ones(size(g_value,1),size(g_value,2));
  if ~isempty(Index_ActivConstr) then
    IneqConstrMatr = [];
    IneqConstrMatr = [df_value' 1; dg_value(Index_ActivConstr,:) Theta(Index_ActivConstr)']; // minus because we have greater than sign in the inequality constraints
  else
    IneqConstrMatr = [df_value' 1]; // minus because we have greater than sign in the inequality constraints
  end
  IneqBoundMatr = [];
  IneqBoundMatr = zeros(length(Index_ActivConstr)+1,1);
  upper_bound = [upper' %inf]';
  lower_bound = [lower' 0]';
  // Optimization by linear programming
  direction = ones(length(x)+1,1); // We add slack variables to the original vector of variables
  
  [direction, lagr, f_opt] = linpro(ObjFuncMatr, IneqConstrMatr, IneqBoundMatr, lower_bound, upper_bound, 0, direction);
    
  // Line search: backtracking algorithm

  // Initialization of t
  //t = - ysk' * d / norm(d)^2;
  t = 0.5; // 2 avant
  if (Log) then
    printf('optim_feasdir: initial value of t = %f\n', t);
  end

  // While the step doesn't verify the armijo condition, we reduce the step
  fk  = _fdir_f(x);
  fsk = _fdir_df(x);
  c1  = 1e-4;
  cg_beta = 0.8;
  Index   = 1;

  while ((Index<ls_ItMX)&((_fdir_f(x+t*direction(1:$-1)) > fk + c1*t*fsk'*direction(1:$-1)) | or(_fdir_g(x+t*direction(1:$-1))>0)))
    t = t * cg_beta;
    Index = Index + 1;
    if (Log) then
      printf('%3d %10.2e %10.2e %10.2e %10.2e\n', Index, t, x - t/cg_beta * direction(1:$-1), _fdir_f(x + t*direction(1:$-1)), d'*_fdir_fs(x + t*direction(1:$-1)));
    end

    if (x_history_defined) then
      x_history($)($+1) = x+t*direction(1:$-1);
    end
  end
  
  x = x + direction(1:$-1)*t;

  if (t<1e-4) then
    if (Log) then
      printf('optim_feasdir: stop on step length: t = %f\n', t);
    end
    x_opt = x;
    return
  end
  
  // If the x step is too small, we stop
  if (norm(x - x_k_m_1)<1.1*XTol) then 
    XTol_Count = XTol_Count + 1;
    if (norm(x - x_k_m_1)>0) then
      XTol = min([XTol norm(x - x_k_m_1)]);
    end
    if Log then
      printf('norm(x - xk) = %f XTol = %f XTol_Count = %d\n', norm(x - x_k_m_1), XTol, XTol_Count);
    end
  else
    XTol_Count = 0;
  end
  if (XTol_Count==MaxMinStepCount) then 
    if (Log) then
      printf('optim_feasdir: stop on tolerance count\n');
    end
    break;
  end
end
x_opt = x;
endfunction

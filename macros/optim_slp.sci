function [x_opt, x_history, ml_history] = optim_slp(slp_f, slp_df, slp_g, slp_dg, slp_h, slp_dh, x0, ...
                                                    ItMX, MaxEvalFunc, delta_ml, upper, lower, Weight, MaxMinStepCount, Log, XTOL, ITOL, ETOL)

if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('MaxMinStepCount','local') then
  MaxMinStepCount = 10;
end

if ~isdef('slp_g','local') | ~isdef('slp_dg','local') then
  NoIneqConstr = %T;
else
  NoIneqConstr = (slp_g==[]);
end

if ~isdef('slp_h','local') | ~isdef('slp_dh','local') then
  NoEqConstr = %T;
else
  NoEqConstr = (slp_h==[]);
end

if ~isdef('slp_f','local') then
  error('optim_slp: slp_f is mandatory');
else
  if typeof(slp_f)=='list' then
    deff('y=_slp_f(x)','y=slp_f(1)(x, slp_f(2:$))');
  else
    deff('y=_slp_f(x)','y=slp_f(x)');
  end
end
if ~isdef('slp_df','local') then
  error('optim_slp: slp_df is mandatory');
else
  if typeof(slp_df)=='list' then
    deff('y=_slp_df(x)','y=slp_df(1)(x, slp_df(2:$))');
  else
    deff('y=_slp_df(x)','y=slp_df(x)');
  end
end

if ~NoEqConstr then
  if typeof(slp_h)=='list' then
    deff('y=_slp_h(x)','y=slp_h(1)(x, slp_h(2:$))');
  else
    deff('y=_slp_h(x)','y=slp_h(x)');
  end

  if typeof(slp_dh)=='list' then
    deff('y=_slp_dh(x)','y=slp_dh(1)(x, slp_dh(2:$))');
  else
    deff('y=_slp_dh(x)','y=slp_dh(x)');
  end
end

if ~NoIneqConstr then
  if typeof(slp_g)=='list' then
    deff('y=_slp_g(x)','y=slp_g(1)(x, slp_g(2:$))');
  else
    deff('y=_slp_g(x)','y=slp_g(x)');
  end

  if typeof(slp_dg)=='list' then
    deff('y=_slp_dg(x)','y=slp_dg(1)(x, slp_vdg(2:$))');
  else
    deff('y=_slp_dg(x)','y=slp_dg(x)');
  end
end

if ~NoIneqConstr then
  g_value  = _slp_g(x0);
  dg_value = _slp_dg(x0);
else
  g_value  = [];
  dg_value = [];
end

if ~NoEqConstr then
  h_value  = _slp_h(x0);
  dh_value = _slp_dh(x0);
else
  h_value  = [];
  dh_value = [];
end

NbIneqConstr = max(size(g_value));
NbEqConstr   = max(size(h_value));

[nargout, nargin] = argn();

if ~isdef('Weight','local') |isempty(Weight) then
  Weight = 1;
end

if ~isdef('XTOL','local') then
  XTOL = 1e-4;
end

if ~isdef('ITOL','local') then
  ITOL = 1e-4;
end

if ~isdef('ETOL','local') then
  ETOL = 1e-4;
end

if length(Weight)~=NbIneqConstr + NbEqConstr then
  Weight = ones(NbIneqConstr + NbEqConstr,1)*Weight(1);
  Weight = Weight / norm(Weight);
end

WeightIneq = [];
WeightEq   = [];

if ~NoIneqConstr then
  WeightIneq = Weight(1:NbIneqConstr);
end
if ~NoEqConstr then
  WeightEq = Weight(NbIneqConstr+1:NbIneqConstr+NbEqConstr);
end

if (length(delta_ml)==1) then
  DeltaMLMax = 20;
  delta_ml_max = ones(size(x0,1), size(x0,2))*DeltaMLMax*delta_ml;
  delta_ml     = ones(size(x0,1), size(x0,2))*delta_ml;
  delta_ml_min = delta_ml;
else
  DeltaMLMax = 20;
  delta_ml_max = DeltaMLMax*delta_ml
  delta_ml_min = delta_ml;
end
x_history_defined = (nargout>=2);
if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end

ml_history_defined = (nargout==3);

// dg must return a matrix NbVar x NbConstr
// x0 must be a column vector of size NbVar
// g must be a column vector of size NbConstr
// df must be a column vector of size NbVar
// delta_ml can either be a scalar or a column vector of size NbVar
// lower and upper must be column vectors of size NbVar

NbConstr = NbEqConstr + NbIneqConstr;

EqConstrMatr = [];
EqBoundMatr  = [];
IneqConstrMatr = [];
IneqBoundMatr  = [];

f_value  = _slp_f(x0);
df_value = _slp_df(x0);
NbVar    = max(size(x0)); 

x = x0;
Reduce_coeff   = 0.3;
RndReduceCoeff = 0.1*Reduce_coeff;
Increase_coeff   = 1.3;
RndIncreaseCoeff = 0.1*Increase_coeff;
MaxML          = 2;
x_k_m_1 = x0;
x_k_m_2 = x0;
x       = x0;

Iteration  = 0;
EvalFunc   = 0;
XTol_Count = 0;

offset = 0;

Debug = %F;

// Verification of the size of all the parameters
if (size(x0,1)~=NbVar) then
  error('x0 must be a column vector of size nbvar');
end
if ~NoIneqConstr then
  if size(g_value,1)~=NbIneqConstr then
    error('g must return a column vector of size nbconstr');
  end
  if ((size(dg_value,1)~=NbVar)&(size(dg_value,2)~=NbIneqConstr)) then
    error('dg must return a matrix of size nbvar x nbconstr');
  end
end
if ~NoEqConstr then
  if size(h_value,1)~=NbEqConstr then
    error('h must return a column vector of size nbconstr');
  end
  if ((size(dh_value,1)~=NbVar)&(size(dh_value,2)~=NbEqConstr)) then
    error('dh must return a matrix of size nbvar x nbconstr');
  end
end
if (size(df_value,1)~=NbVar) then
  error('df must return a column vector of size nbvar');
end
if isempty(delta_ml) then
  // Automatic computation of the Move limits
  // From:
  // K. Y. Chan, S. J. Skerlos, P. Papalambros, "An Adaptative Sequential Linear Programming Algorithm for Optimal Design Problems With Probabilistic Constraints"
  // Trans. Of The ASME, Vol. 129, February 2007
  
  for i=1:NbVar
    delta_ml(i) = 0;
    if ~NoIneqConstr then
      for j=1:NbIneqConstr
        delta_ml(i) = max([delta_ml(i) abs(g_value(j)^2) / sum(abs(dg_value(:,j)).^2)]);
      end
    end
    if ~NoEqConstr then
      for j=1:NbEqConstr
        delta_ml(i) = max([delta_ml(i) abs(h_value(j)^2) / sum(abs(dh_value(:,j)).^2)]);
      end
    end
  end
elseif ((size(delta_ml)~=[1 1]) | (size(delta_ml,1)~=NbVar)) then
  error('delta_ml must be a scalar or a column vector of size nbvar');
end
if (size(lower,1)~=NbVar) then
  error('lower must be a column vector of size NbVar');
end
if (size(upper,1)~=NbVar) then
  error('upper must be a column vector of size NbVar');
end

if or(g_value>0) & Log then
  warning('Initial solution not feasible. Finding a new starting point.');
  printf('constraint value :'); disp(g_value');
  printf('current initial point :'); disp(x');
end

if (ml_history_defined) then
  ml_history = list();
  ml_history($+1) = delta_ml;
end

Iteration = 0;

while((Iteration<ItMX)&(EvalFunc<MaxEvalFunc))
  Iteration = Iteration + 1;
  EvalFunc  = EvalFunc + 1;
  x_k_m_2 = x_k_m_1;
  x_k_m_1 = x;
  // Building local linear models
  // For the objective function
  if (Iteration~=1) then
    f_value  = _slp_f(x);
    df_value = _slp_df(x);
  end
  // For the constraints
  if (Iteration~=1) then
    if ~NoIneqConstr then
      g_value  = _slp_g(x);
      dg_value = _slp_dg(x);
    end
    if ~NoEqConstr then
      h_value  = _slp_h(x);
      dh_value = _slp_dh(x);
    end
  end
  // Initialisation of the slack variables for equality constraints
  w_eq_plus  = zeros(size(h_value,1), size(h_value,2));
  Index = find(h_value>=0);
  w_eq_plus(Index) = abs(h_value(Index));
  w_eq_plus = w_eq_plus;
  w_eq_minus = zeros(size(h_value,1), size(h_value,2));
  Index = find(h_value<0);
  w_eq_minus(Index) = abs(h_value(Index));
  w_eq_minus = w_eq_minus;
  // Initialisation of the slack variables for inequality constraints
  w_ineq = abs(g_value);
  w_ineq = w_ineq;

  // Building the matrix constraint
  // Preparation of the slack variables for equality and inequality constraints
  if ~NoEqConstr & Iteration==1 then
    w_mat_eq_plus  = - diag(WeightEq)*eye(NbEqConstr, NbEqConstr);
    w_mat_eq_minus = diag(WeightEq)*eye(NbEqConstr, NbEqConstr);
  end
  if ~NoIneqConstr & Iteration==1 then
    w_mat_ineq = - diag(WeightIneq)*eye(NbIneqConstr, NbIneqConstr);
  end
  // The Inequality constraints
  if ~NoIneqConstr then
    IneqConstrMatr = [dg_value' w_mat_ineq zeros(NbIneqConstr, NbEqConstr) zeros(NbIneqConstr, NbEqConstr)];
    for i=1:NbIneqConstr
      IneqBoundMatr(i) = - g_value(i) + dg_value(:,i)'*x + offset - w_ineq(i);
    end
  else
    w_mat_ineq     = [];
    IneqConstrMatr = [];
    IneqBoundMatr  = [];
  end
  // The equality constraints
  if ~NoEqConstr then
    EqConstrMatr = [dh_value' zeros(NbEqConstr, NbIneqConstr) w_mat_eq_plus w_mat_eq_minus];
    for i=1:NbEqConstr
      EqBoundMatr(i) = - h_value(i) + dh_value(:,i)'*x + offset - w_eq_plus(i) + w_eq_minus(i);
    end
  else
    EqConstrMatr = [];
    EqBoundMatr  = [];
    w_mat_eq_plus  = [];
    w_mat_eq_minus = [];
  end
  // Optimization by linear programming
  x_extended = [x' w_ineq' w_eq_plus' w_eq_minus']'; // We add slack variables to the original vector of variables
  // First part of the lower / upper bounds are related to move limits and
  // second part are related to slack variables for inequality constraints
  // Sometimes, x+delta is below lower !! This formula prevents this phenomenon
  upper_bound(1:length(x0)) = max([min([x+delta_ml upper],'c') lower], 'c');
  upper_bound(length(x0)+1:length(x_extended)) = %inf; 
//  coeff_w = 1.01;
//  upper_bound(length(x0)+1:length(x0)+length(w_ineq)) = coeff_w*w_ineq; 
//  upper_bound(length(x0)+length(w_ineq)+1:length(x0)+length(w_ineq)+length(w_eq_plus)) = coeff_w*w_eq_plus; 
//  upper_bound(length(x0)+length(w_ineq)+length(w_eq_plus)+1:length(x0)+length(w_ineq)+length(w_eq_plus)+length(w_eq_minus)) = coeff_w*w_eq_minus; 
  // Sometimes, x-delta is above upper !! This formula prevents this phenomenon
  lower_bound(1:length(x0)) = min([max([x-delta_ml lower],'c') upper], 'c');
  lower_bound(length(x0)+1:length(x_extended)) = 0;
  df_value_extended = [df_value' ones(1,length(w_ineq)) ones(1,length(w_eq_plus)) ones(1,length(w_eq_minus))]';

  if (Debug)
    printf('DEBUG: %d / %d\n', Iteration, ItMX);
    if ~NoIneqConstr then
      printf('inequality constraints\n');
      disp(_slp_g(x_extended(1:length(x0))));
      printf('w_ineq\n');
      disp(w_ineq')
    end
    if ~NoEqConstr then
      printf('equality constraints\n');
      disp(_slp_h(x_extended(1:length(x0))));
      printf('w_eq_plus\n');
      disp(w_eq_plus');
      printf('w_eq_minus\n');
      disp(w_eq_minus');
    end
    printf('DEBUGDEBUGDEBUG\n');
    disp(w_eq_plus');
    disp(w_eq_minus');
    disp(w_ineq');
    printf('df_value_extended ');
    disp(size(df_value_extended))
    disp(df_value_extended)
    printf('EqConstrMat ');
    disp(size(EqConstrMatr))
    disp(EqConstrMatr)
    printf('IneqConstrMat ');
    disp(size(IneqConstrMatr))
    disp(IneqConstrMatr)
    printf('EqBoundMatr ');
    disp(size(EqBoundMatr))
    disp(EqBoundMatr)
    printf('IneqboundMatr ');
    disp(size(IneqBoundMatr))
    disp(IneqBoundMatr)
    printf('lower_bound ');
    disp(size(lower_bound))
     printf('upper_bound ');
    disp(size(upper_bound))
    disp(upper_bound - lower_bound)
    printf('x_extended ');
    disp(size(x_extended))
    disp(x_extended)
    printf('Nb Eq Constr = %d Nb Ineq Constr = %d\n', NbEqConstr, NbIneqConstr);
  end

  [x_extended, lagr, f_opt] = linpro(df_value_extended, [EqConstrMatr' IneqConstrMatr']', [EqBoundMatr' IneqBoundMatr']', ...
                                     lower_bound, ...
                                     upper_bound, ...
                                     NbEqConstr, x_extended);
                                     //NbEqConstr, 'v');

//  if (x_extended(1)>upper_bound(1)+%eps) | (x_extended(1)<lower_bound(1)-%eps) | (x_extended(2)>upper_bound(2)+%eps) | (x_extended(2)<lower_bound(2)-%eps)
//    printf('DEBUG: %d / %d\n', Iteration, ItMX);
//    disp(upper_bound')
//    disp(lower_bound')
//    disp(x_extended')
//    disp(df_value_extended)
//    disp([EqConstrMatr' IneqConstrMatr']')
//    disp([EqBoundMatr' IneqBoundMatr']')
//    disp(NbEqConstr)
//    disp(x_extended)
//  end
  
  
  if (Debug) then
    printf('x_extended after\n');
    disp(x_extended)
  end
  x      = x_extended(1:length(x0));
  w_ineq = x_extended(length(x0)+1:length(x0)+NbIneqConstr);
  w_eq_plus  = x_extended(NbIneqConstr+1:NbIneqConstr+NbEqConstr);
  w_eq_minus = x_extended(NbIneqConstr+NbEqConstr+1:NbIneqConstr+NbEqConstr+NbEqConstr);
  
  if (x_history_defined) then
    x_history($+1) = x;
  end
  
  // Update move limits
  if (Debug) then
    printf('delta_ml before'); disp(delta_ml')
  end
  for i=1:length(x)
    if ((x(i) - x_k_m_1(i))*(x_k_m_1(i) - x_k_m_2(i))<0) then
      delta_ml(i) = delta_ml(i) * (Reduce_coeff + (2*rand(1,1)-1)*RndReduceCoeff);
      delta_ml(i) = min([delta_ml(i) delta_ml_max(i)]);
      delta_ml(i) = max([delta_ml(i) delta_ml_min(i)]); // To prevent the collapse of the move limits
    elseif ((x(i) - x_k_m_1(i))*(x_k_m_1(i) - x_k_m_2(i))>0) then
      delta_ml(i) = delta_ml(i) * (Increase_coeff + (2*rand(1,1)-1)*RndIncreaseCoeff);
      delta_ml(i) = min([delta_ml(i) delta_ml_max(i)]);
      delta_ml(i) = max([delta_ml(i) delta_ml_min(i)]); // To prevent the collapse of the move limits
    else
      // No changes in delta_ml
    end
    if (delta_ml(i)>(upper(i) - lower(i)) / (2*MaxML)) then
      delta_ml(i) = (upper(i) - lower(i)) / (2*MaxML);
    end
  end
  if (Debug) then
    printf('delta_ml after'); disp(delta_ml')
  end

  if (ml_history_defined) then
    ml_history($+1) = delta_ml;
  end

  // If the x step is too small, we stop
  if (norm(x - x_k_m_1)<XTOL) then 
    XTol_Count = XTol_Count + 1;
    if Log then
      printf('norm(x - xk) = %f XTol = %f XTol_Count = %d\n', norm(x - x_k_m_1), XTOL, XTol_Count);
    end    
  else
    XTol_Count = 0;
  end

  if ~NoIneqConstr then
    if (max(g_value(find(g_value>0)))<ITOL) then
      x_opt = x;
      return;
    end
  end

  if ~NoEqConstr then
    if (max(abs(h_value))<ETOL) then
      x_opt = x;
      return;
    end
  end

  if (XTol_Count==MaxMinStepCount) then 
    if (Log) then
      printf('optim_slp: stop on tolerance count\n');
    end
    x_opt = x;
    return;
  end
end

x_opt = x;
endfunction

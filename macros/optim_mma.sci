function [x_opt, x_history] = optim_mma(mma_f, mma_df, mma_g, mma_dg, mma_h, mma_dh, x0, upper, lower, ItMX, MaxEvalFunc, XTol, Log)

// Option for compilation of the deff's
option = 'n'; // 'c' for compiled deff. 'n' is better for debugging

[nargout, nargin] = argn();

// Verification of the parameters
if ~isdef('x0','local') then
  error('optim_mma: error,  x0 is mandatory');
end
if ~isdef('upper','local') then
  upper = %inf*ones(size(x0,1), size(x0,2));
end
if ~isdef('lower','local') then
  lower = -%inf*ones(size(x0,1), size(x0,2));
end
if ~isdef('ItMX','local') then
  ItMX = 10;
end
if ~isdef('MaxEvalFunc','local') then
  MaxEvalFunc = ItMX*1000;
end
if ~isdef('XTol','local') then
  XTol = 1e-6;
end
if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('mma_f','local') then
  error('optim_mma: mma_f is mandatory');
else
  if typeof(mma_f)=='list' then
    deff('y=_mma_f(x)','y=mma_f(1)(x, mma_f(2:$))', opt=option);
  else
    deff('y=_mma_f(x)','y=mma_f(x)', opt=option);
  end
end
if ~isdef('mma_df','local') then
  error('optim_mma: mma_df is mandatory');
else
  if typeof(mma_df)=='list' then
    deff('y=_mma_df(x)','y=mma_df(1)(x, mma_df(2:$))', opt=option);
  else
    deff('y=_mma_df(x)','y=mma_df(x)', opt=option);
  end
end
if ~isdef('mma_g','local') then
  g_defined = %F;
else
  if typeof(mma_g)=='list' then
    deff('y=_mma_g(x)','y=mma_g(1)(x, mma_g(2:$))', opt=option);
  else
    deff('y=_mma_g(x)','y=mma_g(x)', opt=option);
  end
  g_defined = %T;
end
if ~isdef('mma_dg','local') then
else
  if typeof(mma_dg)=='list' then
    deff('y=_mma_dg(x)','y=mma_dg(1)(x, mma_dg(2:$))', opt=option);
  else
    deff('y=_mma_dg(x)','y=mma_dg(x)', opt=option);
  end
end
if ~isdef('mma_h','local') then
  h_defined = %F;
else
  if typeof(mma_h)=='list' then
    deff('y=_mma_h(x)','y=mma_h(1)(x, mma_h(2:$))', opt=option);
  else
    deff('y=_mma_h(x)','y=mma_h(x)', opt=option);
  end
  h_defined = %T;
end
if ~isdef('mma_dg','local') then
else
  if typeof(mma_dh)=='list' then
    deff('y=_mma_dh(x)','y=mma_dh(1)(x, mma_dh(2:$))', opt=option);
  else
    deff('y=_mma_dh(x)','y=mma_dh(x)', opt=option);
  end
end

x_history_defined = (nargout==2);
if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end

// x0 must be a column vector
// g must return a column vector of constraint objective function values
Iteration = 0;
EvalFunc  = 0;

NbVars = length(x0);

// Initialisation of the constraints and of the lagrange coefficients
// For inequality constraints
if g_defined then 
  g_value = _mma_g(x0);
  nb_ineq_constr = length(g_value);
  lambda_opt     = ones(1,nb_ineq_constr);
  lambda_lower   = zeros(1,nb_ineq_constr);
  lambda_upper   = %inf*ones(1,nb_ineq_constr);
end
// For equality constraints
if h_defined then 
  h_value = _mma_h(x0);
  nb_eq_constr = length(h_value);
  lambda_opt   = [ones(1,nb_eq_constr)       lambda_opt];
  lambda_lower = [-%inf*ones(1,nb_eq_constr) lambda_lower];
  lambda_upper = [%inf*ones(1,nb_eq_constr)  lambda_upper];
end
//// For bound constraints
//lambda_opt   = [ones(1,NbVars)      lambda_opt];
//lambda_lower = [zeros(1,NbVars)     lambda_lower];
//lambda_upper = [%inf*ones(1,NbVars) lambda_upper];

randomized = %F;

s_mma = 0.7; // for the move limits

_mma_p = zeros(NbVars, nb_ineq_constr + nb_eq_constr); _mma_p_0 = zeros(NbVars,1);
_mma_q = zeros(NbVars, nb_ineq_constr + nb_eq_constr); _mma_q_0 = zeros(NbVars,1);
_mma_r = zeros(NbVars, nb_ineq_constr + nb_eq_constr); _mma_r_0 = zeros(NbVars,1);

x_opt        = x0;
_mma_x_k_m_1 = x0;
_mma_x_k_m_2 = x0;

_mma_u = x_opt + (upper - lower);
_mma_l = x_opt - (upper - lower);

// The move limits
_mma_alpha = 0.9*_mma_l + 0.1*x_opt;
_mma_beta  = 0.1*_mma_u + 0.1*x_opt;
_mma_eta   = s_mma*ones(size(upper,1),size(upper,2));

while ((Iteration<ItMX)&(EvalFunc<MaxEvalFunc))
  _mma_x_k_m_2 = _mma_x_k_m_1;
  _mma_x_k_m_1 = x_opt;
  Iteration = Iteration + 1;
  EvalFunc  = EvalFunc  + 1;
  // Build the MMA model
  // The objective function
  f_value  = _mma_f(x_opt);
  df_value = _mma_df(x_opt);
  _mma_r_0 = f_value;
  for i=1:NbVars
    _mma_p_0(i) = (_mma_u(i) - x_opt(i))^2*(max([0 df_value(i)]));
    _mma_q_0(i) = - (x_opt(i) - _mma_l(i))^2*(min([0 df_value(i)]));
    _mma_r_0    = _mma_r_0 - _mma_p_0(i) / (_mma_u(i) - x_opt(i)) - _mma_q_0(i) / (x_opt(i) - _mma_l(i));
  end
    
  // The inequality constraints
  if g_defined then 
    g_value  = _mma_g(x_opt);
    dg_value = _mma_dg(x_opt);
    for jj=1:nb_ineq_constr
      _mma_r(jj) = g_value(jj);
      // u : NbVars x NbConstraints
      // y : NbConstraints
      // p : NbVars x NbConstraints
      // x : NbVars
      // q : NbVars x NbConstraints
      // l : NbVars x NbConstraints
      for ii=1:NbVars
        _mma_p(ii,jj) = (_mma_u(ii) - x_opt(ii))^2*(max([0 dg_value(ii,jj)]));
        _mma_q(ii,jj) = - (x_opt(ii) - _mma_l(ii))^2*(min([0 dg_value(ii,jj)]));
        _mma_r(jj)   = _mma_r(jj) - _mma_p(ii,jj) / (_mma_u(ii) - x_opt(ii)) - _mma_q(ii,jj) / (x_opt(ii) - _mma_l(ii));
      end
    end
    // First the MMAified constraints (7 first lines)
    // Then the move limits (last lines)
  end

  // The equality constraints
  if h_defined then 
    h_value  = _mma_h(x_opt);
    dh_value = _mma_dh(x_opt);
    for jj=1:nb_eq_constr
      _mma_r(length(g_value)+jj) = h_value(jj);
      // u : NbVars x NbConstraints
      // y : NbConstraints
      // p : NbVars x NbConstraints
      // x : NbVars
      // q : NbVars x NbConstraints
      // l : NbVars x NbConstraints
      for ii=1:NbVars
        _mma_p(ii,length(g_value)+jj) = (_mma_u(ii) - _optx(ii))^2*(max([0 dh_value(ii,jj)]));
        _mma_q(ii,length(g_value)+jj) = - (_optx(ii) - _mma_l(ii))^2*(min([0 dh_value(ii,jj)]));
        _mma_r(length(g_value)+jj)   = _mma_r(length(g_value)+jj) - _mma_p(ii,length(g_value)+jj) / (_mma_u(ii) - _optx(ij)) - _mma_q(ii,length(g_value)+jj) / (_optx(ii) - _mma_l(ii));
      end
    end
    // First the MMAified constraints (7 first lines)
    // Then the move limits (last lines)
  end

  if ~isdef('_to_x','local') then 
    deff('y = _to_x(lambda)','for ii=1:NbVars ...
                                _Aux_1(ii) = _mma_p_0(ii); ...
                                _Aux_2(ii) = _mma_q_0(ii); ...
                                if g_defined then ...
                                  for jj=1:length(g_value) ...
                                    _Aux_1(ii) = _Aux_1(ii) + lambda(jj)*_mma_p(ii,jj); ...
                                    _Aux_2(ii) = _Aux_2(ii) + lambda(jj)*_mma_q(ii,jj); ...
                                  end ...
                                end ...
                                if h_defined then ...
                                  for jj=1:length(h_value) ...
                                    _Aux_1(ii) = _Aux_1(ii) + lambda(length(g_value)+jj)*_mma_p(ii,length(g_value)+jj); ...
                                    _Aux_2(ii) = _Aux_2(ii) + lambda(length(g_value)+jj)*_mma_q(ii,length(g_value)+jj); ...
                                  end ...
                                end ...
                                _Aux_1(ii) = sqrt(_Aux_1(ii)); ...
                                _Aux_2(ii) = sqrt(_Aux_2(ii)); ...
                                y(ii) = (_Aux_1(ii)*_mma_l(ii) + _Aux_2(ii)*_mma_u(ii)) / (_Aux_1(ii) + _Aux_2(ii)); ...
                              end; ...
                              Index = find(y>_mma_u);     y(Index) = _mma_u(Index); ...
                              Index = find(y>_mma_beta);  y(Index) = _mma_beta(Index); ...
                              Index = find(y<_mma_l);     y(Index) = _mma_l(Index); ...
                              Index = find(y<_mma_alpha); y(Index) = _mma_alpha(Index);');
  end
                              
  if ~isdef('f_mma','local') then 
    deff('y = f_mma(lambda)','_x = _to_x(lambda); ...
                              y = _mma_r_0; ...
                              for ii=1:NbVars ...
                                y = y + (_mma_p_0(ii) / (_mma_u(ii) - _x(ii)) + _mma_q_0(ii) / (_x(ii) - _mma_l(ii))); ...
                              end ...
                              if g_defined then ...
                                for jj=1:length(g_value) ...
                                  y = y + lambda(jj)*_mma_r(jj); ...
                                  for ii=1:NbVars ...
                                    y = y + lambda(jj)*(_mma_p(ii,jj) / (_mma_u(ii) - _x(ii)) + _mma_q(ii,jj) / (_x(ii) - _mma_l(ii))); ...
                                  end ...
                                end ...
                              end ...
                              if h_defined then ...
                                for jj=1:length(h_value) ...
                                  y = y + lambda(length(g_value)+jj)*_mma_r(length(g_value)+jj); ...
                                  for ii=1:NbVars ...
                                    y = y + lambda(length(g_value)+jj)*(_mma_p(ii,length(g_value)+jj) / (_mma_u(ii) - _x(ii)) + _mma_q(ii,length(g_value)+jj) / (_x(ii) - _mma_l(ii))); ...
                                  end ...
                                end ...
                              end ...
                              y = -y;');
  end
                       
  if ~isdef('df_mma','local') then 
    deff('y = df_mma(lambda)','y = derivative(f_mma,lambda'')'';');
  end
                          
  if ~isdef('f_mma_optim','local') then
    deff('[f, df, ind_out] = f_mma_optim(lambda, ind_in)','f  = f_mma(lambda); ...
                                                           df = df_mma(lambda); ...
                                                           ind_out = ind_in');
  end
  // Find a solution

  if (Log) then
    printf('optim_mma: Iteration %d / %d\n', Iteration , ItMX);
  end

  lambda_opt = ones(size(lambda_opt,1),size(lambda_opt,2));
  
  [f_opt, lambda_opt] = optim(f_mma_optim, 'b', lambda_lower, lambda_upper, lambda_opt, algo = 'gc', iter = MaxEvalFunc / ItMX);

  // Reconstruction of x_opt
  x_opt = _to_x(lambda_opt);
  
  if (Log) then
    printf('objective function value = %f\n', _mma_f(x_opt));
    if g_defined then
      printf('inequality constraints :'); disp(_mma_g(x_opt));
    end
    if h_defined then
      printf('equality constraints :'); disp(_mma_h(x_opt));
    end
  end
  
  if (x_history_defined) then
    x_history($+1) = x_opt;
  end
  
  // Update rule
  if (Iteration==1 | Iteration==2) then
    _mma_l = x_opt - (upper - lower);
    _mma_u = x_opt + (upper - lower);
  else
    for i=1:NbVars
      if (x_opt(i) - _mma_x_k_m_1(i))*(_mma_x_k_m_1(i) - _mma_x_k_m_2(i))<0 then
        if (Log) then printf('optim_mma: sign change\n'); end
        _mma_eta(i) = s_mma;
      elseif (x_opt(i) - _mma_x_k_m_1(i))*(_mma_x_k_m_1(i) - _mma_x_k_m_2(i))>0 then
        if (Log) then printf('optim_mma: same sign\n'); end
        _mma_eta(i) = 1/s_mma;
      elseif (x_opt(i) - _mma_x_k_m_1(i))*(_mma_x_k_m_1(i) - _mma_x_k_m_2(i))==0 then
        if (Log) then printf('optim_mma: no change\n'); end
        _mma_eta(i) = 1;
      end
    end
    
    if (~randomized) then
      _mma_l = x_opt - _mma_eta.*(_mma_x_k_m_1 - _mma_l);
      _mma_u = x_opt + _mma_eta.*(_mma_u - _mma_x_k_m_1);
    else
      _mma_l = x_opt - ((_mma_eta - _mma_eta/2)*rand(1,1) + _mma_eta/2).*(_mma_x_k_m_1 - _mma_l);
      _mma_u = x_opt + ((_mma_eta - _mma_eta/2)*rand(1,1) + _mma_eta/2).*(_mma_u - _mma_x_k_m_1);
    end
    
    if (x_opt==_mma_x_k_m_1)&(_mma_x_k_m_1==_mma_x_k_m_2) then
      return;
    end
    
    _mma_alpha = 0.9*_mma_l + 0.1*x_opt;
    _mma_beta  = 0.9*_mma_u + 0.1*x_opt;
  end
  if norm(x_opt-_mma_x_k_m_1)<XTol then
    printf('optim_mma: break on XTol\n');
    break;
  end
end
endfunction

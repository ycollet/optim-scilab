function [x_opt, x_history] = optim_LFOPc(lfopc_df, lfopc_g, lfopc_dg, lfopc_h, lfopc_dh, x0, ...
                                          delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start)

// Option for compilation of the deff's
option = 'n'; // 'c' for compiled deff. 'n' is better for debugging

[nargout, nargin] = argn();

x_history_defined = (nargout==2);

if x_history_defined then
  x_history = list();
  x_history($+1) = x0;
end

// Verification of some parameters
if ~isdef('delta_step','local') then
  delta_step = 1;
end
if ~isdef('m','local') then
  m = 5;
end
if ~isdef('delta_inc','local') then
  delta_inc = 0.001;
end
if ~isdef('GradTOL','local') then
  GradTOL = 1e-5;
end
if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('p_start','local') then
  p_start  = 1.01;
end

if ~isdef('lfopc_df','local') then
  error('optim_LFOPc: lfopc_df is mandatory');
else
  if typeof(lfopc_df)=='list' then
    deff('y=_lfopc_df(x)','y=lfopc_df(1)(x, lfopc_df(2:$))', opt=option);
  else
    deff('y=_lfopc_df(x)','y=lfopc_df(x)', opt=option);
  end
end
if ~isdef('x0','local') then
  error('optim_LFOPc: x0 is mandatory');
end

if ~isdef('delta_t','local') then
  // We use a formula to size automatically the value of delta_t
  delta_t = sqrt(delta_step / (5*norm(_lfopc_df(x0))));
end

if ~isdef('lfopc_g','local') then
  g_defined = %F;
  g      = [];
  grad_g = [];
else
  g_defined = (lfopc_g~=[]);
  if (g_defined & typeof(lfopc_g)=='list') then
    deff('y=_lfopc_g(x)','y=lfopc_g(1)(x, lfopc_g(2:$))', opt=option);
  else
    deff('y=_lfopc_g(x)','y=lfopc_g(x)', opt=option);
  end
  if (lfopc_dg~=[] & typeof(lfopc_dg)=='list') then
    deff('y=_lfopc_dg(x)','y=lfopc_dg(1)(x, lfopc_dg(2:$))', opt=option);
  else
    deff('y=_lfopc_dg(x)','y=lfopc_dg(x)', opt=option);
  end
end

if ~isdef('lfopc_h','local') then
  h_defined = %F;
  h      = [];
  grad_h = [];
else
  h_defined = (lfopc_h~=[]);
  if (h_defined & typeof(lfopc_h)=='list') then
    deff('y=_lfopc_h(x)','y=lfopc_h(1)(x, lfopc_h(2:$))', opt=option);
  else
    deff('y=_lfopc_h(x)','y=lfopc_h(x)', opt=option);
  end
  if (lfopc_dh~=[] & typeof(lfopc_dh)=='list') then
    deff('y=_lfopc_dh(x)','y=lfopc_dh(1)(x, lfopc_dh(2:$))', opt=option);
  else
    deff('y=_lfopc_dh(x)','y=lfopc_dh(x)', opt=option);
  end
end

grad_f0 = _lfopc_df(x0);

if (g_defined) then
  g0      = _lfopc_g(x0);
  grad_g0 = _lfopc_dg(x0);
else
  g0      = [];
  grad_g0 = [];
end

if (h_defined) then
  h0      = _lfopc_h(x0);
  grad_h0 = _lfopc_dh(x0);
else
  h0      = [];
  grad_h0 = [];
end

epsilon_h = 1e-3*ones(size(h0,1), size(h0,2));
epsilon_g = 1e-3*ones(size(g0,1), size(g0,2));
epsilon_x = 1e-3;

mu     = 1.0;
_eta   = 1.0;
_beta  = 1.0;
_alpha = 1.0;

NbVars            = length(x0);
NbIneqConstraints = length(g0);
NbEqConstraints   = length(h0);

if (g_defined & ~h_defined) then
  clear lfopc_func;
  deff('y=lfopc_func(x)','grad_f0 = _lfopc_df(x); ...
                          grad_g0 = _lfopc_dg(x); ...
                          g0      = _lfopc_g(x); ...
                          Aux_ineq = abs(grad_g0).^2; Aux_ineq(find(Aux_ineq==0)) = %eps; ...
                          for i=1:NbIneqConstraints ...
                            _alpha(i) = 0; ...
                            for j=1:NbVars ...
                              _alpha(i) = _alpha(i) + (mu *(abs(_eta*grad_f0(j)) + 1) / Aux_ineq(j,i)) * (g0(i)>0); ...
                            end ...
                          end ...
                          Aux_ineq = 0; ...
                          for i=1:NbIneqConstraints ...
                            Aux_ineq = Aux_ineq + _alpha(i) * g0(i) * grad_g0(:,i); ...
                          end ...
                          y = _eta*grad_f0 + 2*Aux_ineq;', opt=option);
end

if (~g_defined & h_defined) then
  clear lfopc_func;
  deff('y=lfopc_func(x)','grad_f0 = _lfopc_df(x); ...
                          grad_h0 = _lfopc_dh(x); ...
                          h0      = _lfopc_h(x); ...
                          Aux_eq  = abs(grad_h0).^2; Aux_eq(find(Aux_eq==0)) = %eps ; ...
                          for i=1:NbEqConstraints ...
                            _beta(i) = 0; ...
                            for j=1:NbVars ...
                              _beta(i) = _beta(i) + mu *(abs(_eta*grad_f0(j)) + 1) / Aux_eq(j,i) ; ...
                            end ...
                          end ...
                          Aux_eq = 0; ...
                          for i=1:NbEqConstraints ...
                            Aux_eq = Aux_eq + _beta(i) * h0(i) * grad_h0(:,i); ...
                          end ...
                          y = _eta*grad_f0 + 2*Aux_eq', opt=option);
end
if (g_defined & h_defined) then
  clear lfopc_func;
  deff('y=lfopc_func(x)','grad_f0 = _lfopc_df(x); ...
                          grad_g0 = _lfopc_dg(x); ...
                          g0      = _lfopc_g(x); ...
                          grad_h0 = _lfopc_dh(x); ...
                          h0      = _lfopc_h(x); ...
                          Aux_ineq = abs(grad_g0).^2; Aux_ineq(find(Aux_ineq==0)) = %eps ; ...
                          for i=1:NbIneqConstraints ...
                            _alpha(i) = 0; ...
                            for j=1:NbVars ...
                              _alpha(i) = _alpha(i) + (mu *(norm(_eta*grad_f0(j)) + 1) / Aux_ineq(j,i)) * (g0(i)>0) ; ...
                            end ...
                          end ...
                          Aux_ineq = 0; ...
                          for i=1:NbIneqConstraints ...
                            Aux_ineq = Aux_ineq + _alpha(i) * g0(i) * grad_g0(:,i); ...
                          end ...
                          Aux_eq = abs(grad_h0).^2; Aux_eq(find(Aux_eq==0)) = %eps ; ...
                          for i=1:NbEqConstraints ...
                            _beta(i) = 0; ...
                            for j=1:NbVars ...
                              _beta(i) = _beta(i) + mu *(abs(_eta*grad_f0(j)) + 1) / Aux_eq(j,i) ; ...
                            end ...
                          end ...
                          Aux_eq = 0; ...
                          for i=1:NbEqConstraints ...
                            Aux_eq = Aux_eq + _beta(i) * h0(i) * grad_h0(:,i); ...
                          end ...
                          y = _eta*grad_f0 + 2*Aux_ineq + 2*Aux_eq', opt=option);
end

// Phase 0
mu      = 1e2;
_eta    = 1;
delta_t = sqrt(delta_step / (5*norm(_lfopc_df(x0))));

if (x_history_defined) then
  [x_tmp, x_history_aux] = optim_LFOP(lfopc_func, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
  if (length(x_history_aux)>2) then
    x_history = lstcat(x_history,x_history_aux);
  else
    if Log then
      printf('optim_LFOPc: Phase 0 terminated without iterations\n');
    end
    x_tmp = x0;
  end
else
  x_tmp = optim_LFOP(lfopc_func, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
end

if (x_tmp~=x0) then
  x_opt = x_tmp;
end

if (g_defined) then 
  g0 = _lfopc_g(x_opt); 
  NbConstr_g = sum(g0<=epsilon_g);
  OverShoot  = sum((g0>epsilon_g) .* g0);
else
  g0         = [];
  NbConstr_g = [];
  OverShoot  = [];
end

if (h_defined) then 
  h0 = _lfopc_h(x_opt); 
  NbConstr_h = sum(h0<=epsilon_h & h0>=-epsilon_h);
else
  h0         = []; 
  NbConstr_h = [];
end

if ((NbConstr_h+NbConstr_g~=0) | OverShoot <= 1e2*epsilon_x) then
  Phase = 1;
else
  Phase = 2;
end

if (Phase==1) then
  // Phase 1
  mu      = 1e4;
  _eta    = 1;
  delta_t = sqrt(delta_step / (5*norm(_lfopc_df(x_opt))));

  if (x_history_defined) then
    [x_tmp, x_history_aux] = optim_LFOP(lfopc_func, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
    if (length(x_history_aux)>2) then
      x_history = lstcat(x_history,x_history_aux);
    else
      if Log then
        printf('optim_LFOPc: Phase 1 terminated without iterations\n');
      end
      x_tmp = x0;
    end
  else
    x_tmp = optim_LFOP(lfopc_func, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
  end
  if (x_tmp~=x0) then
    x_opt = x_tmp;
  end
  Phase = 2;
end

if (Phase==2) then
  // Phase 2
  _eta    = 0;
  delta_t = sqrt(delta_step / (5*norm(_lfopc_df(x_opt))));

  if (g_defined) then
    g0 = _lfopc_g(x_opt);
    NbConstr_g = sum(g0<=epsilon_g);
    OverShoot  = sum((g0>epsilon_g) .* g0);
  else
    g0         = [];
    NbConstr_g = [];
    OverShoot  = [];
  end

  if (h_defined) then
    h0 = _lfopc_h(x_opt);
    NbConstr_h = sum(h0<=epsilon_h & h0>=-epsilon_h);
  else
    h0         = [];
    NbConstr_h = [];
  end
  
  if (x_history_defined) then
    [x_tmp, x_history_aux] = optim_LFOP(lfopc_func, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
    if (length(x_history_aux)>2) then
      x_history = lstcat(x_history,x_history_aux);
    else 
      if Log then
        printf('optim_LFOPc: Phase 2 terminated without iterations\n');
      end
      x_tmp = x0;
    end
  else
    x_tmp = optim_LFOP(lfopc_func, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
  end
  
  if (x_tmp~=x0) then
    x_opt = x_tmp;
  end
end
endfunction

function [x_opt, x_history] = optim_sumt(sumt_f, sumt_df, sumt_g, sumt_dg, sumt_h, sumt_dh, x0, ItMX, NbLoop, Log, r_p, eta, r_max, Algorithm)

[nargout, nargin] = argn();

x_history_defined = (nargout==2);

epsilon = 1e-4;

if ~isdef('sumt_g','local') | ~isdef('sumt_dg','local') then
  g_defined = %F;
  g  = [];
  dg = [];
else
  g_defined = (sumt_g~=[]);
end

if ~isdef('sumt_h','local') | ~isdef('sumt_dh','local') then
  h_defined = %F;
  h  = [];
  dh = [];
else
  h_defined = (sumt_h~=[]);
end

if ~isdef('sumt_f','local') then
  error('optim_sumt: sumt_f is mandatory');
else
  if typeof(sumt_f)=='list' then
    deff('y=_sumt_f(x)','y=sumt_f(1)(x, sumt_f(2:$))');
  else
    deff('y=_sumt_f(x)','y=sumt_f(x)');
  end
end

if ~isdef('sumt_df','local') then
  error('optim_sumt: sumt_df is mandatory');
else
  if typeof(sumt_df)=='list' then
    deff('y=_sumt_df(x)','y=sumt_df(1)(x, sumt_df(2:$))');
  else
    deff('y=_sumt_df(x)','y=sumt_df(x)');
  end
end

if (g_defined) then
  if typeof(sumt_g)=='list' then
    deff('y=_sumt_g(x)','y=sumt_g(1)(x, sumt_g(2:$))');
  else
    deff('y=_sumt_g(x)','y=sumt_g(x)');
  end

  if typeof(sumt_dg)=='list' then
    deff('y=_sumt_dg(x)','y=sumt_dg(1)(x, sumt_dg(2:$))');
  else
    deff('y=_sumt_dg(x)','y=sumt_dg(x)');
  end
end

if (h_defined) then
  if typeof(sumt_h)=='list' then
    deff('y=_sumt_h(x)','y=sumt_h(1)(x, sumt_h(2:$))');
  else
    deff('y=_sumt_h(x)','y=sumt_h(x)');
  end

  if typeof(sumt_dh)=='list' then
    deff('y=_sumt_dh(x)','y=sumt_dh(1)(x, sumt_dh(2:$))');
  else
    deff('y=_sumt_dh(x)','y=sumt_dh(x)');
  end
end

lambda = 0.1;

// Initialization of variables
f_value  = _sumt_f(x0);
df_value = _sumt_df(x0);
if (g_defined) then
  g_value  = _sumt_g(x0);
  dg_value = _sumt_dg(x0);
  lambda_g = lambda*ones(size(g_value,1),size(g_value,2));
end
if (h_defined) then
  h_value  = _sumt_h(x0);
  dh_value = _sumt_dh(x0);
  lambda_h = lambda*ones(size(h_value,1),size(h_value,2));
end

// Normalisation par le gradient

if g_defined then
  df_value_aux = norm(df_value);
  for i=1:size(dg_value,2)
    norm_g(i,1) = df_value_aux ./ max([norm(dg_value(:,i)) %eps]);
  end
end

if h_defined then
  df_value_aux = norm(df_value);
  for i=1:size(dh_value,2)
    norm_h(i,1) = df_value_aux ./ max([norm(dh_value(:,i)) %eps]);
  end
end

if ~isdef('Log','local') then
  Log = %F;
end
if isempty(Log) then
  Log = %F;
end

if ~isdef('NbLoop','local') then
  NbLoop = 3;
end
if isempty(NbLoop) then
  NbLoop = 3;
end

if ~isdef('r_p','local') then
  // Set penalty parameters:
  // Augmented Lagrangian function is of the form:
  // AL = F(X) + \SUM_i P(LAMBDA_i,C_i(X),RHO),
  // where
  // P(LAMBDA,Y,RHO) = Y ( LAMBDA + 0.5 * RHO * Y ),  if C_i(X) is an equality constraint or LAMBDA + RHO * Y > 0, and
  // P(LAMBDA,Y,RHO) = - 0.5 * LAMBDA^2 / RHO, otherwise.
  // Assuming that LAMBDA_i = 0 for all i, it is clear that  P(LAMBDA_i,C_i(X),RHO) = 0.5 * RHO *C_i(X)^2
  // And that the valud of THO that balances F(X) and \SUM_i P(LAMBDA_i,C_i(X),RHO) is given by
  // RHO = F(X) / ( 0.5 * \SUM_i C_i(X)^2 ),
  // Where the sum is made over the violated constraints
  
  aux = 0;
  if g_defined then
    aux = aux + sum(g_value(find(g_value>0)).^2);
  end
  if h_defined then
    aux = aux + sum(h_value.^2);
  end
  if aux==0 then
    r_p = 10.0;
  else
    r_p = 10 * max(1, abs(f_value)) / (0.5*aux);
    r_p = max(1e-6,min(10.0, r_p));
  end
end
if isempty(r_p) then
  aux = 0;
  if g_defined then
    aux = aux + sum(g_value(find(g_value>0)).^2);
  end
  if h_defined then
    aux = aux + sum(h_value.^2);
  end
  if aux==0 then
    r_p = 10.0;
  else
    r_p = 10 * max(1, abs(f_value)) / (0.5*aux);
    r_p = max(1e-6,min(10.0, r_p));
  end
end

if ~isdef('r_max','local') then
  r_max = 100*r_p;
end
if isempty(r_max) then
  r_max = 100*r_p;
end

if ~isdef('eta','local') then
  eta = abs(r_max / r_p)^(1/NbLoop);
end
if isempty(eta) then
  eta = abs(r_max / r_p)^(1/NbLoop);
end

if ~isdef('Algorithm','local') then
  Algorithm = 'gc';
end
if isempty(Algorithm) then
  Algorithm = 'gc';
end

if Log |%T then
  printf('optim_sumt: r_p = %.4f / r_max = %.4f / eta = %.4f / NbLoop = %d\n', r_p, r_max, eta, NbLoop); 
end

global EvalFunc

EvalFunc = 0;

if (g_defined & ~h_defined) then
  // A penalized function for inequality constraints
  deff('[f_lap, df_lap, ind] = obj_func(x, ind)','global EvalFunc; ...
                                                  EvalFunc = EvalFunc + 1; ...
                                                  g_value = norm_g .* _sumt_g(x); dg_value = _sumt_dg(x); ...
                                                  f_value = _sumt_f(x); df_value = _sumt_df(x); ...
                                                  for i=1:length(g_value) ...
                                                    g_value_prim(i) = max([g_value(i) -lambda_g(i)/(2*r_p)]); ...
                                                  end; ...
                                                  f_lap = f_value + sum(lambda_g.*g_value_prim) + r_p*sum(g_value_prim.^2); ...
                                                  df_lap = df_value; ...
                                                  for i=1:size(dg_value,1) ...
                                                    Aux(i) = 0; ...
                                                    for j=1:size(dg_value,2) ...
                                                      if g_value(j)>(-lambda_g(j)/(2*r_p)) then ...
                                                        Aux(i) = Aux(i) + lambda_g(j) * norm_g(j) * dg_value(i,j) + r_p * 2 * norm_g(j) * dg_value(i,j) * g_value(j); ...
                                                      end; ...
                                                    end ...
                                                  end; ...
                                                  df_lap = df_lap + Aux;');
end

if (~g_defined & h_defined) then
  // A penalized function for equality constraints
  deff('[f_lap, df_lap, ind] = obj_func(x, ind)','df_lap = 0; ...
                                                  global EvalFunc; ...
                                                  EvalFunc = EvalFunc + 1; ...
                                                  f_value = _sumt_f(x); df_value = _sumt_df(x); ...
                                                  h_value = norm_h .* _sumt_h(x); dh_value = _sumt_dh(x); ...
                                                  f_lap = f_value + sum(lambda_h.*h_value) + r_p*sum(h_value.^2); ...
                                                  df_lap = df_value; ...
                                                  for i=1:size(dh_value,1) ...
                                                    Aux(i) = 0; ...
                                                    for j=1:size(dh_value,2) ...
                                                      Aux(i) = Aux(i) + lambda_h(j) * norm_h(j) * dh_value(i,j) + r_p * 2 * norm_h(j) * dh_value(i,j) * h_value(j); ...
                                                    end; ...
                                                  end; ...
                                                  df_lap = df_lap + Aux;');
end

if (g_defined & h_defined) then
  // A penalized function for inequality and equality constraints
  deff('[f_lap, df_lap, ind] = obj_func(x, ind)','df_lap = 0; ...
                                                  global EvalFunc; ...
                                                  EvalFunc = EvalFunc + 1; ...
                                                  g_value = norm_g .* _sumt_g(x); dg_value = _sumt_dg(x); ...
                                                  f_value = _sumt_f(x); df_value = _sumt_df(x); ...
                                                  h_value = norm_h .* _sumt_h(x); dh_value = _sumt_dh(x); ...
                                                  for i=1:length(g_value) ...
                                                    g_value_prim(i) = max([g_value(i) -lambda_g(i)/2*r_p]); ...
                                                  end; ...
                                                  f_lap = f_value + sum(lambda_g.*g_value_prim) + r_p*sum(g_value_prim.^2); ...
                                                  df_lap = df_value; ...
                                                  for i=1:size(dg_value,1) ...
                                                    Aux(i) = 0; ...
                                                    for j=1:size(dg_value,2) ...
                                                      if g_value(j)>(-lambda_g(j)/(2*r_p)) then ...
                                                        Aux(i) = Aux(i) + lambda_g(j) * norm_g(j) * dg_value(i,j) + r_p * 2 * norm_g(j) * dg_value(i,j) * g_value(j); ...
                                                      end; ...
                                                    end ...
                                                  end; ...
                                                  df_lap = df_lap + Aux; ...
                                                  Aux = []; ...
                                                  f_lap = f_lap + sum(lambda_h.*h_value) + r_p*sum(h_value.^2); ...
                                                  for i=1:size(dh_value,1) ...
                                                    Aux(i) = 0; ...
                                                    for j=1:size(dh_value,2) ...
                                                      Aux(i) = Aux(i) + lambda_h(j) * norm_h(j) * dh_value(i,j) + r_p * 2 * norm_h(j) * dh_value(i,j) * h_value(j); ...
                                                    end; ...
                                                  end; ...
                                                  df_lap = df_lap + Aux;');
end

x_opt = x0;

if x_history_defined then
  x_history = list();
  x_history($+1) = x_opt;
end

for i=1:NbLoop
  [f_opt, x_opt] = optim(obj_func, x_opt, algo=Algorithm, 'ar', ItMX / NbLoop, ItMX / NbLoop, 1e-6, 1e-6);
  if (g_defined) then
    g_value = norm_g .* _sumt_g(x_opt);
    for j=1:length(lambda_g)
      lambda_g(j) = lambda_g(j) + 2*r_p*max([g_value(j), -lambda_g(j)/(2*r_p)]);
    end
  end
  if (h_defined) then
    h_value = norm_h .* _sumt_h(x_opt);
    for j=1:length(lambda_h)
      lambda_h(j) = lambda_h(j) + 2*r_p*h_value(j);
    end
  end
  
  r_p = r_p * eta;
  
  if (r_p>r_max) then r_p = r_max; end
  
  if (x_history_defined) then
    x_history($+1) = x_opt;
  end
  if (Log) then
    printf('Iteration %d / %d - r_p = %f, f_obj = %f:\n', i, NbLoop, r_p, sumt_f(x_opt));
    if g_defined then
      printf('max overshoot for constraints (must be all negative or null): %f\n', max(sumt_g(x_opt)));
    end
    if h_defined then
      printf('max overshoot for constraints (must be all null): %f\n', max(sumt_h(x_opt)));
    end
  end
end
endfunction

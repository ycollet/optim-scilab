function [x_opt, x_history, ml_history] = optim_dynq(dynq_f, dynq_df, dynq_g, dynq_dg, dynq_h, dynq_dh, x0, x_low, x_up, NbOuterIter, GradTOL, XTol, MaxEvalFunc, Log, lambda_init)
// Parameter verification
[nargout, nargin] = argn();

x_history_defined  = (nargout>=2);
ml_history_defined = (nargout==3);

if x_history_defined then
  x_history = list();
end

if ml_history_defined then
  ml_history      = list();
end

if ~isdef('dynq_g','local') | ~isdef('dynq_dg','local') then
  g_defined = %F;
  g  = [];
  dg = [];
else
  g_defined = (dynq_g~=[]);
end

if ~isdef('dynq_h','local') | ~isdef('dynq_dh','local') then
  h_defined = %F;
  h  = [];
  dh = [];
else
  h_defined = (dynq_h~=[]);
end

if ~isdef('NbOuterIter','local') then
  NbOuterIter = 100;
end
if ~isdef('MaxEvalFunc','local') then
  MaxEvalFunc = 100;
end
if ~isdef('GradTOL','local') then
  GradTOL = 1e-5;
end
if ~isdef('Log','local') then
  Log = %F;
end

if ~isdef('dynq_f','local') then
  error('optim_dynq: dynq_f is mandatory');
else
  if typeof(dynq_f)=='list' then
    deff('y=_dynq_f(x)','y=dynq_f(1)(x, dynq_f(2:$))');
  else
    deff('y=_dynq_f(x)','y=dynq_f(x)');
  end
end

if ~isdef('dynq_df','local') then
  error('optim_dynq: dynq_df is mandatory');
else
  if typeof(dynq_df)=='list' then
    deff('y=_dynq_df(x)','y=dynq_df(1)(x, dynq_df(2:$))');
  else
    deff('y=_dynq_df(x)','y=dynq_df(x)');
  end
end

if (g_defined) then
  if typeof(dynq_g)=='list' then
    deff('y=_dynq_g(x)','y=dynq_g(1)(x, dynq_g(2:$))');
  else
    deff('y=_dynq_g(x)','y=dynq_g(x)');
  end

  if typeof(dynq_dg)=='list' then
    deff('y=_dynq_dg(x)','y=dynq_dg(1)(x, dynq_dg(2:$))');
  else
    deff('y=_dynq_dg(x)','y=dynq_dg(x)');
  end
end

if (h_defined) then
  if typeof(dynq_h)=='list' then
    deff('y=_dynq_h(x)','y=dynq_h(1)(x, dynq_h(2:$))');
  else
    deff('y=_dynq_h(x)','y=dynq_h(x)');
  end

  if typeof(dynq_dh)=='list' then
    deff('y=_dynq_dh(x)','y=dynq_dh(1)(x, dynq_dh(2:$))');
  else
    deff('y=_dynq_dh(x)','y=dynq_dh(x)');
  end
end

if ~isdef('lambda_init','local') then
  lambda_init = [];
end

f_value  = _dynq_f(x0);
df_value = _dynq_df(x0);
if (g_defined) then
  g_value  = _dynq_g(x0);
  dg_value = _dynq_dg(x0);

  for i=1:size(dg_value,2)
    lambda_g(i,1) = norm(df_value) / max([norm(dg_value(:,i)) %eps]);
  end
else
  g_value  = [];
  dg_value = [];
end

if (h_defined) then
  h_value  = _dynq_h(x0);
  dh_value = _dynq_dh(x0);
  for i=1:size(dh_value,2)
    lambda_h(i,1) = norm(df_value) / max([norm(dh_value(:,i)) %eps]);
  end
else
  h_value  = [];
  dh_value = [];
end

//// Inhibit the weighting of constraints
if g_defined then lambda_g = ones(size(lambda_g,1),size(lambda_g,2)); end
if h_defined then lambda_h = ones(size(lambda_h,1),size(lambda_h,2)); end

Algorithm  = 'gc';
InitMatrix = 1;
LambdaInit = 1;
MaxLambda  = %inf;
MinLambda  = 0;

// Boundary limits and starting point
upper_bound = [];
lower_bound = [];
lambda_opt  = [];

// We add the lambdas
// For the inequality constraints
if g_defined then
  lower_bound = MinLambda * ones(length(g_value),1);
  upper_bound = MaxLambda * ones(length(g_value),1);
  lambda_opt  = LambdaInit * ones(length(g_value),1);
end
// For the equality constraints
if h_defined then
  lower_bound = [lower_bound' -MaxLambda * ones(1,length(h_value))]';
  upper_bound = [upper_bound'  MaxLambda * ones(1,length(h_value))]';
  lambda_opt  = [lambda_opt'  LambdaInit * ones(length(h_value))]';
end
// For the low bound constraint
lower_bound = [lower_bound' MinLambda * ones(1,length(x0))]';
upper_bound = [upper_bound' MaxLambda * ones(1,length(x0))]'
lambda_opt  = [lambda_opt' LambdaInit * ones(1,length(x0))]';
// For the up bound constraint
lower_bound = [lower_bound' MinLambda * ones(1,length(x0))]';
upper_bound = [upper_bound' MaxLambda * ones(1,length(x0))]'
lambda_opt  = [lambda_opt' LambdaInit * ones(1,length(x0))]';

if ~isempty(lambda_init) then
  lambda_opt = lambda_init;
end

mode_ieee = ieee();
ieee(2);

x_opt   = x0;
x_k_m_1 = x0;
x_k_m_2 = x0;

ItMX = ceil(MaxEvalFunc / NbOuterIter);

if g_defined then n_ineq_constr = length(g_value); end
if h_defined then n_eq_constr   = length(h_value); end

// Parameters for move limits
eta       = ones(size(x_opt,1),size(x_opt,2));
delta_max = (x_up - x_low) / 2.0;
delta_min = (x_up - x_low) / 100.0;
delta     = (delta_max - delta_min) / 4 + delta_min;
eta_inc   = 2;
eta_dec   = 1 / eta_inc;
Index_eta_inc = 0;
Index_eta_dec = 0;
Max_eta_inc = 2;
Max_eta_dec = 2;

Debug = %F;

if size(lambda_opt,2)~=1 then lambda_opt = lambda_opt'; end

for It=1:NbOuterIter
  x_k_m_2 = x_k_m_1;
  x_k_m_1 = x_opt;
  // Building the quadratic model
  // For the objective function
  f_value_old  = f_value;
  f_value      = _dynq_f(x_k_m_1);
  df_value_old = df_value;
  df_value     = _dynq_df(x_k_m_1);
  
  if (g_defined) then
    g_value_old  = g_value;
    g_value      = _dynq_g(x_k_m_1);
    dg_value_old = dg_value;
    dg_value     = _dynq_dg(x_k_m_1);

    for j=1:size(dg_value,2)
      lambda_g(j,1) = norm(df_value) / max([norm(dg_value(:,j)) %eps]);
    end
  end

  if (h_defined) then
    h_value_old  = h_value;
    h_value      = _dynq_h(x_k_m_1);
    dh_value_old = dh_value;
    dh_value     = _dynq_dh(x_k_m_1);

    for j=1:size(dh_value,2)
      lambda_h(j,1) = norm(df_value) / max([norm(dh_value(:,j)) %eps]);
    end
  end

  // Inhibit the weighting of constraints
  if g_defined then lambda_g = ones(size(lambda_g,1),size(lambda_g,2)); end
  if h_defined then lambda_h = ones(size(lambda_h,1),size(lambda_h,2)); end

  x_delta = x_k_m_2 - x_k_m_1;
  Index_delta = find(x_delta<=10*%eps & x_delta>=-10*%eps);
  for i=1:length(Index_delta)
    if sign(x_delta(i))==0 then
      x_delta(i) = 10*%eps;
    else
      x_delta(i) = sign(x_delta(i))*10*%eps;
    end
  end

  if (It==1) then
    a = InitMatrix*ones(length(x0),1);
  else
    a = ((df_value_old - df_value) ./ (x_delta));
//    Index_tmp = find(isnan(a)==%T);
//    a(Index_tmp) = 1;
  end
  
  // For the constraints
  if g_defined then
    for j=1:n_ineq_constr
      if (It==1) then
        b(:,j) = InitMatrix*ones(length(x0),1);
      else
        b(:,j) = ((dg_value_old(:,j) - dg_value(:,j)) ./ (x_delta));
      end
//      Index_tmp = find(isnan(b(:,j))==%T);
//      b(Index_tmp,j) = 1;
    end
  end

  if h_defined then
    for j=1:n_eq_constr
      if (It==1) then
        c(:,j) = InitMatrix*ones(length(x0),1);
      else
        c(:,j) = ((dh_value_old(:,j) - dh_value(:,j)) ./ (x_delta));
      end
//      Index_tmp = find(isnan(c(:,j))==%T);
//      c(Index_tmp,j) = 1;
    end
  end
  
  if Log & Debug then  
    printf('a,b,c coefficients = ');
    disp(a);
    if g_defined then disp(b'); end
    if h_defined then disp(c'); end
  end
  //
  // Extraction of the optimal value of x given a value of lambda
  //

  //
  // A faire
  // Vérifier l'insertion de d_g (voir dans le document lyx)
  // reprendre la methode de la page 365 de MMA pour l'affectation des bornes.
  
  // dans lambda, on a:
  // lambda(1:n_ineq_constr)                     : les coeff de lagrange des contraintes d'inegalite
  // lambda(n_ineq_constr+1:n_ineq_constr+n_dim) : les coeff de lagrange des contraintes de borne inf
  // lambda(n_ineq_constr+n_dim+1:$)             : les coeff de lagrange des contraintes de borne sup
  
  deff('_x = to_x(lambda)','for ii=1:length(x_opt) ...
                              _num(ii) = df_value(ii) - a(ii)*x_k_m_1(ii);...
                              _den(ii) = a(ii);...
                              if g_defined then ...
                                for jj=1:length(g_value) ...
                                  _num(ii) = _num(ii) + lambda(jj)*lambda_g(jj)*(dg_value(ii,jj) - b(ii,jj)*x_k_m_1(ii));...
                                  _den(ii) = _den(ii) + lambda(jj)*lambda_g(jj)*b(ii,jj);...
                                end;...
                              end;...
                              if h_defined then ...
                                for jj=1:length(h_value) ...
                                  _num(ii) = _num(ii) + lambda(length(g_value)+jj)*lambda_h(jj)*(dh_value(ii,jj) - c(ii,jj)*x_k_m_1(ii));...
                                  _den(ii) = _den(ii) + lambda(length(g_value)+jj)*lambda_h(jj)*c(ii,jj);...
                                end;...
                              end;...
                              _num(ii) = _num(ii) - 2*lambda(length(g_value)+length(h_value)+ii)*x_k_m_1(ii);...
                              _den(ii) = _den(ii) + 2*lambda(length(g_value)+length(h_value)+ii);...
                            end; ...
                            _x = - _num ./ _den; ...
                            for ii=1:length(x_opt) ...
                              if _den(ii)<=0 then ...
                                _Aux_1(ii) = f_value + df_value(ii)*(x_low_ml(ii) - x_k_m_1(ii)) + 0.5*a(ii)*(x_low_ml(ii) - x_k_m_1(ii))^2;...
                                _Aux_2(ii) = f_value + df_value(ii)*(x_up_ml(ii)  - x_k_m_1(ii)) + 0.5*a(ii)*(x_up_ml(ii)  - x_k_m_1(ii))^2;...
                                if g_defined then ...
                                  for jj=1:length(g_value) ...
                                    _Aux_1(ii) = _Aux_1(ii) + lambda(jj)*lambda_g(jj)*(g_value(jj) + dg_value(ii,jj)*(x_low_ml(ii) - x_k_m_1(ii)) + 0.5*b(ii,jj)*(x_low_ml(ii) - x_k_m_1(ii))^2); ...
                                    _Aux_2(ii) = _Aux_2(ii) + lambda(jj)*lambda_g(jj)*(g_value(jj) + dg_value(ii,jj)*(x_up_ml(ii)  - x_k_m_1(ii)) + 0.5*b(ii,jj)*(x_up_ml(ii)  - x_k_m_1(ii))^2); ...
                                  end;...
                                end;...
                                if h_defined then ...
                                  for jj=1:length(h_value) ...
                                    _Aux_1(ii) = _Aux_1(ii) + lambda(length(g_value)+jj)*lambda_h(jj)*(h_value(jj) + dh_value(ii,jj)*(x_low_ml(ii) - x_k_m_1(ii)) + 0.5*c(ii,jj)*(x_low_ml(ii) - x_k_m_1(ii))^2); ...
                                    _Aux_2(ii) = _Aux_2(ii) + lambda(length(g_value)+jj)*lambda_h(jj)*(h_value(jj) + dh_value(ii,jj)*(x_up_ml(ii)  - x_k_m_1(ii)) + 0.5*c(ii,jj)*(x_up_ml(ii)  - x_k_m_1(ii))^2); ...
                                  end;...
                                end;...
                                _Aux_1(ii) = _Aux_1(ii) + lambda(length(g_value)+length(h_value)+ii)*((x_low_ml(ii) - x_k_m_1(ii))^2 - (((x_up_ml(ii) - x_low_ml(ii))/2)^2)); ...
                                _Aux_2(ii) = _Aux_2(ii) + lambda(length(g_value)+length(h_value)+ii)*((x_up_ml(ii)  - x_k_m_1(ii))^2 - (((x_up_ml(ii) - x_low_ml(ii))/2)^2)); ...
                                if (_Aux_1(ii)>_Aux_2(ii)) then ...
                                  _x(ii) = x_up_ml(ii); ...
                                else ...
                                  _x(ii) = x_low_ml(ii); ...
                                end ...
                              end ...
                            end; ...
                            Index = find(_x>x_up_ml); ...
                            _x(Index) = x_up_ml(Index);...
                            Index = find(_x<x_low_ml); ...
                            _x(Index) = x_low_ml(Index);','n');

  //
  // Expression of the objective function
  //

  deff('y=dynq2_opt_f(lambda)','_x = to_x(lambda);...
                                y = f_value + sum(df_value.*(_x - x_k_m_1)) + 0.5*sum(a.*(_x - x_k_m_1).^2);...
                                if g_defined then ...
                                  for ii=1:length(g_value) ...
                                    y = y + lambda(ii)*lambda_g(ii)*(g_value(ii) + sum(dg_value(:,ii).*(_x - x_k_m_1)) + 0.5*sum(b(:,ii).*(_x - x_k_m_1).^2));...
                                  end ...
                                end ...
                                if h_defined then ...
                                  for ii=1:length(h_value) ...
                                    y = y + lambda(length(g_value)+ii)*lambda_h(ii)*(h_value(ii) + sum(dh_value(:,ii).*(_x - x_k_m_1)) + 0.5*sum(c(:,ii).*(_x - x_k_m_1).^2));...
                                  end ...
                                end ...
                                y = y + sum(lambda(length(g_value)+length(h_value)+1:length(g_value)+length(h_value)+length(x_k_m_1)) .*((_x - x_k_m_1).^2  - ((x_up_ml - x_low_ml)/2).^2)); ...
                                y = - y;');

  //
  // Expression of the derivative of the objective function
  //

  deff('y=dynq2_opt_df(lambda)','y = derivative(dynq2_opt_f, lambda);');

  //
  // Gathering objective function and derivative so as to be compatible with the 'optim' function of scilab
  //

  deff('[y,dy,ind_out]=optim_dynq2_opt_f(lambda,ind_in)','y  = dynq2_opt_f(lambda); ...
                                                          dy = dynq2_opt_df(lambda); ...
                                                          ind_out = ind_in;');
  
  // Solving
  if (Log) then
    printf('optim_dynq2: iteration %d / %d\n', It, NbOuterIter);
  end
    
  x_low_ml = max([x_k_m_1 - delta x_low],'c');
  x_up_ml  = min([x_k_m_1 + delta x_up],'c');

  lambda_opt = LambdaInit * ones(size(lambda_opt,1),size(lambda_opt,2));
  
  [f_opt, lambda_opt] = optim(optim_dynq2_opt_f, 'b', lower_bound, upper_bound, lambda_opt, algo=Algorithm, 'ar', ItMX, ItMX, 1e-16, 1e-16,imp=3);
  
  x_opt_old = x_opt;
  x_opt = to_x(lambda_opt);
  
  if Log & Debug then
    printf('lambda_opt = '); disp(lambda_opt');
    if g_defined then printf('lambda_g = '); disp(lambda_g'); end
    if h_defined then printf('lambda_h = '); disp(lambda_h'); end
  end
  
  x_opt_aux = x_opt;
  x_opt = max([x_opt x_low_ml],'c');
  if or(x_opt~=x_opt_aux) & Log then 
    printf('x_opt low saturation\n');
    printf('x_low_ml     = '); disp(x_low_ml');
    printf('x_opt before = '); disp(x_opt_aux');
    printf('x_opt after  = '); disp(x_opt');
  end
  
  x_opt_aux = x_opt;
  x_opt = min([x_opt x_up_ml],'c');
  if or(x_opt~=x_opt_aux) & Log then
    printf('x_opt high saturation\n');
    printf('x_up_ml      = '); disp(x_up_ml');
    printf('x_opt before = '); disp(x_opt_aux');
    printf('x_opt after  = '); disp(x_opt');
  end

  // Updating the move limits
  for ii=1:length(x_opt)
    if (x_opt(ii) - x_k_m_1(ii))*(x_k_m_1(ii) - x_k_m_2(ii))<0 then
      if Log then printf('optim_dynq2: sign change\n'); end
      if Index_eta_inc>=Max_eta_inc then
        if Log then printf('optim_dynq2: sign change - increase\n'); end
        eta(ii,1) = eta_dec;
        Index_eta_inc = 0;
      else
        Index_eta_inc = Index_eta_inc + 1;
      end
    elseif (x_opt(ii) - x_k_m_1(ii))*(x_k_m_1(ii) - x_k_m_2(ii))>0 then
      if Log then printf('optim_dynq2: same sign\n'); end
      if Index_eta_dec>=Max_eta_dec then
        if Log then printf('optim_dynq2: sign change - decrease\n'); end
        eta(ii,1) = eta_inc;
        Index_eta_dec = 0;
      else
        Index_eta_dec = Index_eta_dec + 1;
      end
    elseif (x_opt(ii) - x_k_m_1(ii))*(x_k_m_1(ii) - x_k_m_2(ii))==0 then
      if Log then printf('optim_dynq2: no change\n'); end
      eta(ii,1) = 1;
    end
  end
  
  delta_old = delta;
  delta = delta .* eta;
  
  if Log & Debug then
    printf('optim_dynq2: eta   ='); disp(eta');
    printf('optim_dynq2: delta ='); disp(delta');
    printf('optim_dynq2: delta_old ='); disp(delta_old');
    printf('optim_dynq2: x_opt ='); disp(x_opt');

  end
  
  Index = find(delta>delta_max);
  delta(Index) = delta_max(Index);
  Index = find(delta<delta_min);
  delta(Index) = delta_min(Index);

  if ml_history_defined then
    ml_history($+1)  = list();
    ml_history($)(1) = delta_old;
    ml_history($)(2) = x_k_m_1;
  end

  if (Log) then
    aux = f_value + sum(df_value.*(x_opt - x_k_m_1)) + 0.5*sum(a.*(x_opt - x_k_m_1).^2);
    printf('optim_dynq2: value of the objective function: %f\n',aux);

    if (g_defined) then
      for ii=1:length(g_value) 
        aux = g_value(ii) + sum(dg_value(:,ii).*(x_opt - x_k_m_1)) + 0.5*sum(b(:,ii).*(x_opt - x_k_m_1).^2);
        printf('optim_dynq2: value of the inequality constraint %d = %f\n',ii,aux);
      end
    end
    if (h_defined) then
      for ii=1:length(h_value) 
        aux = h_value(ii) + sum(dh_value(:,ii).*(x_opt - x_k_m_1)) + 0.5*sum(c(:,ii).*(x_opt - x_k_m_1).^2);
        printf('optim_dynq2: value of the equality constraint %d = %f\n',ii,aux);
      end
    end
    for ii=1:length(x_k_m_1)
      aux = x_opt(ii) - x_up_ml(ii);
      printf('optim_dynq2: value of the up bound constraints %d = %f\n', ii, aux);
    end
    for ii=1:length(x_k_m_1)
      aux = x_low_ml(ii) - x_opt(ii);
      printf('optim_dynq2: value of the low bound constraints %d = %f\n', ii, aux);
    end
  end

  if (x_history_defined) then
    x_history($+1) = x_opt;
  end
  
  if ((norm(x_opt_old - x_opt)<XTol) & It>1) then
    if Log then
      printf('optim_dynq2: X tolerance reached. Exit\n');
      ieee(mode_ieee);
      return;
    end
  end
  if ((norm(x_delta)<XTol) & It>1) then
    if Log then
      printf('optim_dynq2: X_delta tolerance reached. Exit\n');
      ieee(mode_ieee);
      return;
    end
  end
end

ieee(mode_ieee);
endfunction

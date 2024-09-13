function [x_opt, x_history] = optim_sqsd(sqsd_f, sqsd_df, x0, ItMX, XTol, GradTol, rho, Log)
// Parameter verification
[nargout, nargin] = argn();

x_history_defined = (nargout==2);

if x_history_defined then
  x_history = list();
end

if ~isdef('ItMX','local') then
  ItMX = 100;
end
if ~isdef('GradTol','local') then
  GradTol = 1e-6;
end
if ~isdef('XTol','local') then
  XTol = 1e-6;
end
if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('rho','local') then
  rho = 1e-3;
end

if ~isdef('sqsd_f','local') then
  error('optim_sqsd: sqsd_f is mandatory');
else
  if typeof(sqsd_f)=='list' then
    deff('y=_sqsd_f(x)','y=sqsd_f(1)(x, sqsd_f(2:$))');
  else
    deff('y=_sqsd_f(x)','y=sqsd_f(x)');
  end
end

if ~isdef('sqsd_df','local') then
  error('optim_sqsd: sqsd_df is mandatory');
else
  if typeof(sqsd_df)=='list' then
    deff('y=_sqsd_df(x)','y=sqsd_df(1)(x, sqsd_df(2:$))');
  else
    deff('y=_sqsd_df(x)','y=sqsd_df(x)');
  end
end

if ~isdef('x0','local') then
  error('optim_sqsd: x0 is mandatory');
end

f_k  = _sqsd_f(x0);
df_k = _sqsd_df(x0);
f_k_m_1 = f_k;
df_k_m_1 = df_k;

x_k     = x0;
x_k_m_1 = x0;

ck = norm(df_k) / rho;

for i=1:ItMX
  if Log then
    printf('optim_sqsd: Iteration %d / %d - f = %f\n',i , ItMX, f_k);
  end
  
  if x_history_defined then
    x_history($+1) = x_k;
  end
  
  x_k_m_1 = x_k;

  if norm(df_k)<GradTol then
    if Log then
      printf('optim_sqsd: stop on GradTol threshold\n');
    end
    
    x_opt = x_k;
    return;
  else
    x_k = x_k_m_1 - df_k / ck;
  end
  
  if norm(df_k)>rho then
    x_k = x_k - rho * df_k / norm(df_k);
  end
  
  if norm(x_k - x_k_m_1)<XTol then
    if Log then
      printf('optim_sqsd: stop on XTol threshold\n');
    end
    
    x_opt = x_k;
    return;
  end
  
  df_k = _sqsd_df(x_k);
  f_k = _sqsd_f(x_k);
  
  ck = 2*(f_k_m_1 - f_k - df_k'*(x_k_m_1 - x_k))/norm(x_k - x_k_m_1)^2
  
  if ck <= 0 then
    ck = 1e-60;
  end
end

x_opt = x_k;
endfunction

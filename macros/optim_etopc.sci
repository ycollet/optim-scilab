function [x_opt, x_history] = optim_etopc(etopc_f, etopc_df, etopc_g, etopc_dg, etopc_h, etopc_dh, x0, ItMX, ItOut, dt, rho, TOL, GradTOL, XTOL, typecg, Log);

[nargout, nargin] = argn();

// inputs
if ~isdef('dt','local') then
  dt = 0.5;
end
if ~isdef('rho','local') then
  rho = 1000;
end
if ~isdef('typecg','local') then
  typecg  = 'fr'; // Other choice: 'pr'
end
if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('x0','local') | ~isdef('etopc_f','local') | ~isdef('etopc_df','local') then
  error('optim_etopc: error, x0, etopc_f and etopc_df are mandatory');
end
if ~isdef('ItMX','local') then
  ItMX = 1000;
end
if ~isdef('ItOut','local') then
  ItOut = 6;
end
if ~isdef('FTOL','local') then
  FTOL = 1e-6;
end
if ~isdef('GradTOL','local') then
  GradTOL = 1e-6;
end
if ~isdef('XTOL','local') then
  XTOL = 1e-6;
end

ineq_constr_defined = isdef('etopc_g','local') & isdef('etopc_dg','local');
eq_constr_defined   = isdef('etopc_h','local') & isdef('etopc_dh','local');

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history         = list();
  x_history($+1)    = list();
  x_history($)($+1) = x0;
end

if (Log) then
  if typecg=='fr' then
     printf('optim_etopc: using Fletcher-Reeves...\n');
  elseif typecg=='pr' then
     printf('optim_etopc: using Polak-Ribiere...\n');
  else
     error('optim_etoc: specify ''fr'' or ''pr'' for typecg...');
  end
end

if ineq_constr_defined & eq_constr_defined then
  deff('y = df_etopc(x)','y = etopc_df(x); ...
                          ineq_constr  = etopc_g(x); ...
                          dineq_constr = etopc_dg(x); ...
                          eq_constr    = etopc_h(x); ...
                          deq_constr   = etopc_dh(x); ...
                          for j=1:length(x) ...
                            for i=1:length(eq_constr) ...
                              y(j) = y(j) + 2*rho*deq_constr(j,i)*eq_constr(i); ...
                            end ...
                            for i=1:length(ineq_constr) ...
                              Aux = max([0 ineq_constr(i)]); ...
                              if (Aux>0) then ...
                                y(j) = y(j) + 2*rho*dineq_constr(j,i)*Aux; ...
                              end ...
                            end ...
                          end');
elseif ineq_constr_defined & ~eq_constr_defined then
  deff('y = df_etopc(x)','y = etopc_df(x); ...
                          ineq_constr  = etopc_g(x); ...
                          dineq_constr = etopc_dg(x); ...
                          for j=1:length(x) ...
                            for i=1:length(ineq_constr) ...
                              Aux = max([0 ineq_constr(i)]); ...
                              if (Aux>0) then ...
                                y(j) = y(j) + 2*rho*dineq_constr(j,i)*Aux; ...
                              end ...
                            end ...
                          end');
elseif ~ineq_constr_defined & eq_constr_defined then
  deff('y = df_etopc(x)','y = etopc_df(x); ...
                          eq_constr    = etopc_h(x); ...
                          deq_constr   = etopc_dh(x); ...
                          for j=1:length(x) ...
                            for i=1:length(eq_constr) ...
                              y(j) = y(j) + 2*rho*deq_constr(j,i)*eq_constr(i); ...
                            end ...
                          end');
elseif ~ineq_constr_defined & ~eq_constr_defined then
  deff('y = df_etopc(x)','y = etopc_df(x);');
  
  ItOut = 1;
end

if Log then
   printf('optim_etopc: running...\n');
end

x_opt = x0;

for i=1:ItOut
//  dt = 0.5 / rho;
  if Log then
    printf('optim_etopc: extern iteration %d / %d\n', i, ItOut);
  end
  
  [x_opt, x_history_aux] = optim_etop(df_etopc, x_opt, dt, ItMX, typecg, XTOL, GradTOL, Log);

  if (x_history_defined) then
    x_history($+1) = list();
    x_history($)   = lstcat(x_history($), x_history_aux);
  end

  rho = rho*10;
end
endfunction

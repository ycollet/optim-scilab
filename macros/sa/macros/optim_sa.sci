function [x_best, f_best, mean_list, var_list, temp_list, f_history, x_history] = optim_sa(x0, sa_f, ItExt, ItInt, T0, Log, ...
                                                                                           temp_law, param_temp_law, ...
                                                                                           neigh_func, param_neigh_func)
// Simulated annealing
// x0         : initial solution
// f          : objective function
// ItExt      : number of temperature decrease
// ItInt      : number of iterations during one temperature step
// T0         : initial temperature
// Log        : print some message during the run of the optimization
// temp_law   : function which controls the temperature decrease
// param_temp_law   : a variable (any kind of structure) which contains the parameters of the temperature law function
// neigh_func       : function which returns a neighbor of a given point
// param_neigh_func : a variable (any kind of structure) which contains the parameters of the neighborhood function

[nargout, nargin] = argn();

if (~isdef('temp_law','local')) then
  temp_law = temp_law_default;
end
if (~isdef('param_temp_law','local')) then
  param_temp_law = [];
end

if (~isdef('neigh_func','local')) then
  neigh_func = neigh_func_default;
end
if (~isdef('param_neigh_func','local')) then
  param_neigh_func = [];
end

if (~isdef('Log','local')) then
  Log = %F;
end

if (nargout>=6) then
  f_history_defined = %T;
  f_history = list();
else
  f_history_defined = %F;
end

if (nargout>=5) then
  temp_list_defined = %T;
  temp_list = [];
else
  temp_list_defined = %F;
end

if (nargout>=7) then
  x_history_defined = %T;
  x_history = list();
else
  x_history_defined = %F;
end

if ~isdef('sa_f','local') then
  error('optim_sa: sa_f is mandatory');
else
  if typeof(sa_f)=='list' then
    deff('y=_sa_f(x)','y=sa_f(1)(x, sa_f(2:$))');
  else
    deff('y=_sa_f(x)','y=sa_f(x)');
  end
end

T = T0;

// Some variables needed to record the behavior of the SA
var_list  = [];
mean_list = [];
temp_list = [];

x_current = x0;
f_current = _sa_f(x_current);

x_best = x_current;
f_best = f_current;

for i=1:ItExt
  f_list = [];
  x_list = list();
  for j=1:ItInt
    x_neigh = neigh_func(x_current,T,param_neigh_func);
    f_neigh = _sa_f(x_neigh);
    if ((f_neigh<=f_current)|(exp(-(f_neigh-f_current)/T)>rand(1,1))) then
      x_current = x_neigh;
      f_current = f_neigh;
    end
    
    f_list = [f_list f_current];
    
    if (f_best>f_current) then
      x_best = x_current;
      f_best = f_current;
    end
    
    if (x_history_defined) then
      x_list($+1) = x_current;
    end
  end

  if (temp_list_defined) then
    temp_list = [temp_list T];
  end
  if (x_history_defined) then
    x_history($+1) = x_list;
  end
  if (f_history_defined) then
    f_history($+1) = f_list;
  end

  // Computation of step_mean and step_var
  step_mean = mean(f_list);
  step_var  = stdev(f_list);
  mean_list = [mean_list step_mean];
  var_list  = [var_list step_var];
  
  if (Log) then
    printf('optim_sa: Temperature step %d / %d - T = %f, E(f(T)) = %f var(f(T)) = %f f_best = %f\n', i, ItExt, T, step_mean, step_var, f_best);
  end

  T = temp_law(T, step_mean, step_var, i, max(size(x_current)), param_temp_law);
end
endfunction


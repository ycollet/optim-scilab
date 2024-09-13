function [x_best, temp_list, x_history] = optim_st(x0, st_f, ItMax, T, st_proba, Log, neigh_func, param_neigh_func)
// Simulated tempering
// x0         : initial solution
// f          : objective function
// ItMax      : number of iterations
// T          : vector of temperatures
// st_proba   : probability of increasing or decreasing the temperature
// Log        : print some message during the run of the optimization
// neigh_func       : function which returns a neighbor of a given point
// param_neigh_func : a variable (any kind of structure) which contains the parameters of the neighborhood function

[nargout, nargin] = argn();

if (~isdef('neigh_func','local')) then
  neigh_func = neigh_func_default;
end
if (~isdef('param_neigh_func','local')) then
  param_neigh_func = [];
end
if (~isdef('Log','local')) then
  Log = %F;
end
if (~isdef('st_proba','local')) then
  st_proba = 0.05;
end

if (nargout>=2) then
  temp_list_defined = %T;
  temp_list = [];
else
  temp_list_defined = %F;
end

if (nargout==3) then
  x_history_defined = %T;
  x_history = list();
else
  x_history_defined = %F;
end

if ~isdef('st_f','local') then
  error('optim_st: st_f is mandatory');
else
  if typeof(st_f)=='list' then
    deff('y=_st_f(x)','y=st_f(1)(x, st_f(2:$))');
  else
    deff('y=_st_f(x)','y=st_f(x)');
  end
end
// Some variables needed to record the behavior of the SA
temp_list = [];
x_best = [];

Index_T   = length(T);
T_neigh   = T(Index_T);
T_current = T_neigh;
x_current = x0;
f_current = _st_f(x_current);
x_best    = x0;
f_best    = f_current;
f_threshold = f_best;

for i=1:ItMax
  if rand(1,1)<st_proba then
    if (rand(1,1) < 0.5) then
      Index_T_neigh = min([Index_T+1 length(T)]);
    else
      Index_T_neigh = max([Index_T-1 1]);
    end

    T_neigh = T(Index_T_neigh);
    if (exp(-(f_current - f_threshold)*(1/T_current - 1/T_neigh))>rand(1,1)) then
      if Log then
        printf('optim_st: temperature transition - from %f to %f\n', T_current, T_neigh);
      end
      T_current = T_neigh;
      Index_T   = Index_T_neigh;
    end
  else
    x_neigh = neigh_func(x_current,T_current,param_neigh_func);
    f_neigh = _st_f(x_neigh);
    if ((f_neigh<=f_current)|(exp(-(f_neigh-f_current)/T_current)>rand(1,1))) then
      x_current = x_neigh;
      f_current = f_neigh;
    end

    if (f_best>f_current) then
      x_best = x_current;
      f_best = f_current;
      f_threshold = f_best;
    end
  end

  if (temp_list_defined) then
    temp_list = [temp_list T_current];
  end
  if (x_history_defined) then
    x_history($+1) = x_current;
  end
  
  if (Log) then
    printf('optim_st: Iteration %d / %d - T = %f, f_best = %f\n', i, ItMax, T_current, f_best);
  end
end
endfunction


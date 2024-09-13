function [x_best, x_history] = optim_pt(x0, pt_f, ItMax, T, pt_proba, Log, neigh_func, param_neigh_func)
// ParallelS tempering
// x0         : initial solution
// f          : objective function
// ItMax      : number of iterations
// T          : vector of temperatures
// pt_proba   : probability of increasing or decreasing the temperature
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
if (~isdef('pt_proba','local')) then
  pt_proba = 0.05;
end

if (nargout>=2) then
  temp_list_defined = %T;
  temp_list = [];
else
  temp_list_defined = %F;
end

if (nargout==2) then
  x_history_defined = %T;
  x_history = list();
  for i=1:length(T)
    x_history(i) = list();
  end
else
  x_history_defined = %F;
end

if ~isdef('pt_f','local') then
  error('optim_pt: pt_f is mandatory');
else
  if typeof(pt_f)=='list' then
    deff('y=_pt_f(x)','y=pt_f(1)(x, pt_f(2:$))');
  else
    deff('y=_pt_f(x)','y=pt_f(x)');
  end
end

x_current = list();

for i=1:length(T)
  x_current(i) = x0;
  f_current(i) = _pt_f(x_current(i));
end

x_best    = x0;
f_best    = f_current(1);
f_threshold = f_best;

for i=1:ItMax
  for j=1:length(T)
    if rand(1,1)<pt_proba then
      if (rand(1,1) < 0.5) then
        Index_T_neigh = min([j+1 length(T)]);
      else
        Index_T_neigh = max([j-1 1]);
      end

      if (exp(-(f_current(j) - f_threshold)*(1/T(j) - 1/T(Index_T_neigh)))>rand(1,1)) then
        if Log then
          printf('optim_pt: temperature transition - from %f to %f\n', T(j), T(Index_T_neigh));
        end
        x_current(j) = x_current(Index_T_neigh);
        f_current(j) = f_current(Index_T_neigh);
      end
    else
      x_neigh = neigh_func(x_current(j),T(j),param_neigh_func);
      f_neigh = _pt_f(x_neigh);
      if ((f_neigh<=f_current(j))|(exp(-(f_neigh-f_current(j))/T(j))>rand(1,1))) then
        x_current(j) = x_neigh;
        f_current(j) = f_neigh;
      end

      if (f_best>f_current(j)) then
        x_best = x_current(j);
        f_best = f_current(j);
        f_threshold = f_best;
      end
    end

    if (x_history_defined) then
      x_history(j)($+1) = x_current(j);
    end

    if (Log) then
      printf('optim_pt: Iteration %d / %d - f_best = %f\n', i, ItMax, f_best);
    end
  end
end
endfunction


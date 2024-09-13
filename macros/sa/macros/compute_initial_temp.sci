function T_init = compute_initial_temp(x0, f, proba_init, ItMX, neigh_func, param_neigh_func)

if (~isdef('neigh_func','local')) then
  neigh_func = neigh_func_default;
end
if (~isdef('param_neigh_func','local')) then
  param_neigh_func = [];
end

f_list    = [];
x_current = x0;
f_current = f(x_current);
f_list    = [f_list f_current];

for i=1:ItMX
  x_current = neigh_func(x_current, 0, param_neigh_func);
  f_current = f(x_current);
  f_list = [f_list f_current];
end

NbInc = 0;
f_sum = 0;

for i=2:size(f_list,2)
  if (f_list(i-1)<f_list(i)) then
    NbInc = NbInc + 1;
    f_sum = f_sum + f_list(i);
  end
end

if (NbInc>0) then
  f_sum = f_sum / NbInc;
end

// proba_init = exp(-delta_f/T_init) -> -delta_f / log(proba_init) = T_init
T_init = - f_sum ./ log(proba_init);
endfunction

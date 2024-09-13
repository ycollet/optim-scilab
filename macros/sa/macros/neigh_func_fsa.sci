function x_neigh = neigh_func_fsa(x_current, T, param)
if (~isdef('param','local') | isempty(param)) then
  param = ones(size(x_current,1),size(x_current,2));
end
x_neigh = x_current + T*param.*tan(%pi*(rand(size(x_current,1),size(x_current,2)) - 0.5));
endfunction

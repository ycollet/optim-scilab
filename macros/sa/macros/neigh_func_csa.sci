function x_neigh = neigh_func_csa(x_current, T, param)
if (~isdef('param','local') | isempty(param)) then
  param = ones(size(x_current,1), size(x_current,2));
end
x_neigh = x_current + param.*sqrt(2)*T*rand(1,1,'norm');
endfunction

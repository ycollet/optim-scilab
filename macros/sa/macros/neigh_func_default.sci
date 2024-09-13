function x_neigh = neigh_func_default(x_current, T, param)
if (~isdef('param','local') | isempty(param)) then
  sa_min_delta = -0.1*ones(size(x_current,1),size(x_current,2));
  sa_max_delta = 0.1*ones(size(x_current,1),size(x_current,2));
else
  sa_min_delta = param(:,1);
  sa_max_delta = param(:,2);
end
x_neigh = x_current + (sa_max_delta - sa_min_delta).*rand(size(x_current,1),size(x_current,2)) + sa_min_delta;
endfunction

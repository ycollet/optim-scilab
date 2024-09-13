function x_neigh = neigh_func_vfsa(x_current, T, param)
if (~isdef('param','local') | isempty(param)) then
  Delta = 0.1*ones(size(x_current,1), size(x_current,2));
else
  Delta = param;
end
y = zeros(size(x_current,1), size(x_current,2));
for i=1:length(y)
  a = rand(1,1);
  y(i) = sign(a-0.5)*T*((1+1/T)^abs(2*a-1)-1);
  x_neigh(i) = x_current(i) + y(i)*Delta(i);
end
endfunction

function T = temp_law_huang(T, step_mean, step_var, temp_stage, n,param)
if (~isdef('param','local')) then
  param = []; // First create the empty param var
end
if (isempty(param)) then
  lambda = 0.01;
else
  lambda = param(1);
end
T = T * exp(-lambda*T/step_var);
endfunction

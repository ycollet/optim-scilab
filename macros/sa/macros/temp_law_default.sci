function T = temp_law_default(T, step_mean, step_var, temp_stage, n,param)
if (~isdef('param','local')) then
  param = []; // First create the empty var param
end
if (isempty(param)) then
  _alpha = 0.9;
else
  _alpha = param(1);
end
T = _alpha*T;
endfunction

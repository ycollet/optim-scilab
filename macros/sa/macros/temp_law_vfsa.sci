function T = temp_law_vfsa(T, step_mean, step_var, temp_stage, n, param)
if (~isdef('param','local')) then
  param = []; // First create the empty param var
end
if (isempty(param)) then
  c = 0.01;
else
  c = param(1);
end
T = T * exp(-c*((temp_stage+1)^(1/n)-(temp_stage)^(1/n)));
endfunction

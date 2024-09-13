function T = temp_law_fsa(T, step_mean, step_var, temp_stage, n, param)
T = T * (1+temp_stage)/(2+temp_stage);
endfunction

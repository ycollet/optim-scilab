function [Crossed_Indiv1, Crossed_Indiv2] = crossover_func_default(Indiv1,Indiv2,param)
  if ~isdef('param','local') then
    param = [];
  end
  
  // We deal with some parameters to take into account the boundary of the domain and the neighborhood size
  if is_param(param,'beta') then
    Beta = get_param(param,'beta');
  else
    Beta = 0;
  end
  
  if is_param(param,'minbound') then
    MinBounds = get_param(param,'minbound');
  else
    MinBounds = -2*ones(size(Indiv1,1),size(Indiv1,2));
  end

  if is_param(param,'maxbound') then
    MaxBounds = get_param(param,'maxbound');
  else
    MaxBounds = 2*ones(size(Indiv1,1),size(Indiv1,2));
  end

  mix = (1 + 2*Beta)*rand(1,1) - Beta;
  Crossed_Indiv1 =     mix*Indiv1 + (1-mix)*Indiv2;
  Crossed_Indiv2 = (1-mix)*Indiv1 +     mix*Indiv2;
  
  Crossed_Indiv1 = max(min(Crossed_Indiv1, MaxBounds),MinBounds);
  Crossed_Indiv2 = max(min(Crossed_Indiv2, MaxBounds),MinBounds);
endfunction

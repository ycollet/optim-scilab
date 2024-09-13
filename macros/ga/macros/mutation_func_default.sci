function Mut_Indiv = mutation_func_default(Indiv,param)
  if ~isdef('param','local') then
    param = [];
  end
  // We deal with some parameters to take into account the boundary of the domain and the neighborhood size
  if is_param(param,'delta') then
    Delta = get_param(param,'delta');
  else
    Delta = 0.1;
  end
  
  if is_param(param,'minbound') then
    MinBounds = get_param(param,'minbound');
  else
    MinBounds = -2*ones(size(Indiv,1),size(Indiv,2));
  end

  if is_param(param,'maxbound') then
    MaxBounds = get_param(param,'maxbound');
  else
    MaxBounds = 2*ones(size(Indiv,1),size(Indiv,2));
  end
  
  Mut_Indiv = Indiv + 2*Delta*rand(size(Indiv,1),size(Indiv,2)) - Delta*ones(size(Indiv,1),size(Indiv,2));
  
  Mut_Indiv = max(min(Mut_Indiv, MaxBounds),MinBounds);
endfunction

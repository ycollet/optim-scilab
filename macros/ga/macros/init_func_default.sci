function Pop_init = init_func_default(popsize,param)
  if ~isdef('param','local') then
    param = [];
  end
  
  // We deal with some parameters to take into account the boundary of the domain and the neighborhood size
  if is_param(param,'dimension') then
    Dim = get_param(param,'dimension');
  else
    Dim = 2;
  end
  
  if is_param(param,'minbound') then
    MinBounds = get_param(param,'minbound');
  else
    MinBounds = -2*ones(1,Dim);
  end

  if is_param(param,'maxbound') then
    MaxBounds = get_param(param,'maxbound');
  else
    MaxBounds = 2*ones(1,Dim);
  end

  // Pop_init must be a list()
  Pop_init = list();
  for i=1:popsize
    Pop_init(i) = (MaxBounds - MinBounds).*rand(size(MaxBounds,1),size(MaxBounds,2)) + MinBounds;
  end
endfunction

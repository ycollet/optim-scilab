function [pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_ga(ga_f, pop_size, nb_generation, p_mut, p_cross, init_func, param, Log, crossover_func, mutation_func, codage_func, selection_strategy, nb_couples, pressure)

[nargout, nargin] = argn();

if (~isdef('codage_func','local')) then
  codage_func = codage_identity;
end

if ~isdef('ga_f','local') then
  error('optim_ga: ga_f is mandatory');
else
  if typeof(ga_f)=='list' then
    deff('y=_ga_f(x)','y=ga_f(1)(x, ga_f(2:$))');
  else
    deff('y=_ga_f(x)','y=ga_f(x)');
  end
end

if ~isdef('pop_size','local') then
  pop_size = 100;
end
if ~isdef('nb_generation','local') then
  nb_generation = 10;
end
if ~isdef('p_mut','local') then
  p_mut = 0.1;
end
if ~isdef('p_cross','local') then
  p_cross = 0.1;
end
if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('selection_strategy','local') then
  selection_strategy = 'elitist';
end
if ~isdef('nb_couples','local') then
  nb_couples = 100;
end
if ~isdef('pressure','local') then
  pressure = 0.05;
end

// Initialization of the population
if (Log) then
  printf('optim_ga: Initialization of the population\n');
end

Pop = list();
Pop = init_func(pop_size);
// Code the individuals
Pop = codage_func(Pop,'code');

for i=1:length(Pop)
  FObj_Pop(i) = _ga_f(Pop(i));
end

if (nargout==4) then
  pop_init = Pop;
  fobj_pop_init = FObj_Pop;
end

FObj_Pop_Max = max(FObj_Pop);
FObj_Pop_Min = min(FObj_Pop);

// Normalization of the efficiency

Efficiency = (1 - pressure) * (FObj_Pop_Max - FObj_Pop)/max([FObj_Pop_Max - FObj_Pop_Min, %eps]) + pressure;

// The genetic algorithm
for i=1:nb_generation
  if (Log) then
    printf('optim_ga: iteration %d / %d', i, nb_generation);
  end
  //
  // Selection
  //
  Indiv1 = list();
  Indiv2 = list();
  Wheel = cumsum(Efficiency);
  for j=1:nb_couples
    // Selection of the first individual in the couple
    Shoot = rand(1,1)*Wheel($);
    Index = 1;
    while((Wheel(Index)<Shoot)&(Index<length(Wheel))) 
      Index = Index + 1;
    end
    Indiv1(j) = Pop(Index);
    // Selection of the second individual in the couple
    Shoot = rand(1,1)*Wheel($);
    Index = 1;
    while((Wheel(Index)<Shoot)&(Index<length(Wheel))) 
      Index = Index + 1;
    end
    Indiv2(j) = Pop(Index);
  end
  //
  // Crossover
  //
  for j=1:nb_couples
    if (p_cross>rand(1,1)) then
      [x1, x2] = crossover_func(Indiv1(j), Indiv2(j));
      Indiv1(j) = x1;
      Indiv2(j) = x2;
    end
  end
  //
  // Mutation
  //
  for j=1:nb_couples
    if (p_mut>rand(1,1)) then
      x1 = mutation_func(Indiv1(j));
      Indiv1(j) = x1;
    end
    if (p_mut>rand(1,1)) then
      x2 = mutation_func(Indiv2(j));
      Indiv2(j) = x2;
    end
  end
  //
  // Computation of the objective functions
  //
  for j=1:nb_couples
    FObj_Indiv1(j) = _ga_f(Indiv1(j));
    FObj_Indiv2(j) = _ga_f(Indiv2(j));
  end
  //
  // Recombination
  //
  select selection_strategy
  case 'elitist' then
    Total_Pop  = lstcat(Pop, Indiv1, Indiv2);
    Total_FObj = [FObj_Pop' FObj_Indiv1' FObj_Indiv2']';
    // Normalization of the efficiency
    FObj_Pop_Max = max(Total_FObj);
    FObj_Pop_Min = min(Total_FObj);
    
    Efficiency = (1 - pressure) * (FObj_Pop_Max - Total_FObj)/max([FObj_Pop_Max - FObj_Pop_Min, %eps]) + pressure;
    [Efficiency, Index_Sort] = sort(Efficiency);
    Efficiency = Efficiency(1:pop_size);
    // Extraction and selection of the phenotype
    Total_FObj = Total_FObj(Index_Sort);
    FObj_Pop   = Total_FObj(1:pop_size);
    // Extraction and selection of the genotype
    Total_Pop = list(Total_Pop(Index_Sort));
    Pop       = list(Total_Pop(1:pop_size));
    Total_Pop = list(); // Reinitialisation of Total_Pop
    Sorted_Total_Pop = list(); // Reinitialisation of Sorted_Total_Pop
    if (Log) then
      printf(' - min / max value found = %f / %f\n', FObj_Pop_Min, FObj_Pop_Max);
    end
  else
    error('optim_ga: wrong selection strategy');
  end
end

pop_opt  = Pop;
Pop_opt  = codage_func(pop_opt,'decode');
fobj_pop_opt = FObj_Pop;
endfunction

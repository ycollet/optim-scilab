function [pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_moga(ga_f, pop_size, nb_generation, p_mut, p_cross, init_func, param, Log, crossover_func, mutation_func, codage_func, selection_strategy, nb_couples, pressure)

[nargout, nargin] = argn();

if ~isdef('ga_f','local') then
  error('optim_moga: ga_f is mandatory');
end

if ~isdef('pop_size','local') then
  pop_size = 100;
end
if ~isdef('nb_generation','local') then
  nb_generation = 10;
end
if ~isdef('p_mut','local') then
  p_mut = 0.01;
end
if ~isdef('p_cross','local') then
  p_cross = 0.7;
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
  printf('optim_moga: Initialization of the population\n');
end

Total_Pop = list();
Total_Pop = init_func(pop_size);
// Code the individuals
Total_Pop = codage_func(Total_Pop,'code');

for i=1:length(Total_Pop)
  Aux_Total_FObj(i,:) = ga_f(Total_Pop(i));
end

// Compute the domination rank
for i=1:size(Aux_Total_FObj,1)
  Index = 0;
  for j=1:size(Aux_Total_FObj,1)
    Index = Index + double(and(Aux_Total_FObj(i,:)<=Aux_Total_FObj(j,:)) & or(Aux_Total_FObj(i,:)<Aux_Total_FObj(j,:)));
  end      
  Total_FObj(i) = (Index + 1);
end

FObj_Pop_Max = max(Total_FObj);
FObj_Pop_Min = min(Total_FObj);

// Normalization of the efficiency

//Efficiency = (1 - pressure) * (FObj_Pop_Max - Total_FObj) / max([FObj_Pop_Max - FObj_Pop_Min %eps]) + pressure;
Efficiency = (1 - pressure) * (Total_FObj - FObj_Pop_Min) / max([FObj_Pop_Max - FObj_Pop_Min %eps]) + pressure;

if (nargout==4) then
  pop_init = Total_Pop;
  fobj_pop_init = Aux_Total_FObj;
end

// The genetic algorithm
for i=1:nb_generation
  if (Log) then
    printf('optim_moga: iteration %d / %d', i, nb_generation);
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
    Indiv1(j) = Total_Pop(Index);
    // Selection of the second individual in the couple
    Shoot = rand(1,1)*Wheel($);
    Index = 1;
    while((Wheel(Index)<Shoot)&(Index<length(Wheel))) 
      Index = Index + 1;
    end
    Indiv2(j) = Total_Pop(Index);
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
  for j=1:length(Indiv1)
    Aux_FObj_Indiv1(j,:) = ga_f(Indiv1(j)); 
    Aux_FObj_Indiv2(j,:) = ga_f(Indiv2(j));
  end

  // Compute the domination rank
  for j=1:size(Aux_FObj_Indiv1,1)
    // We compute the rank for Indiv1
    Index1 = 0; Index2 = 0; Index3 = 0;
    for k=1:size(Aux_FObj_Indiv1,1)
      Index1 = Index1 + double(and(Aux_FObj_Indiv1(j,:)<=Aux_FObj_Indiv1(k,:)) & or(Aux_FObj_Indiv1(j,:)<Aux_FObj_Indiv1(k,:)));
      Index2 = Index2 + double(and(Aux_FObj_Indiv1(j,:)<=Aux_FObj_Indiv2(k,:)) & or(Aux_FObj_Indiv1(j,:)<Aux_FObj_Indiv2(k,:)));
    end
    for k=1:size(Aux_Total_FObj,1)
      Index3 = Index3 + double(and(Aux_FObj_Indiv1(j,:)<=Aux_Total_FObj(k,:)) & or(Aux_FObj_Indiv1(j,:)<Aux_Total_FObj(k,:)));
    end
    FObj_Indiv1(j) = (Index1 + Index2 + Index3 + 1);

    // We compute the rank for Indiv2
    Index1 = 0; Index2 = 0; Index3 = 0;
    for k=1:size(Aux_FObj_Indiv1,1)
      Index1 = Index1 + double(and(Aux_FObj_Indiv2(j,:)<=Aux_FObj_Indiv1(k,:)) & or(Aux_FObj_Indiv2(j,:)<Aux_FObj_Indiv1(k,:)));
      Index2 = Index2 + double(and(Aux_FObj_Indiv2(j,:)<=Aux_FObj_Indiv2(k,:)) & or(Aux_FObj_Indiv2(j,:)<Aux_FObj_Indiv2(k,:)));
    end
    for k=1:size(Aux_Total_FObj,1)
      Index3 = Index3 + double(and(Aux_FObj_Indiv2(j,:)<=Aux_Total_FObj(k,:)) & or(Aux_FObj_Indiv2(j,:)<Aux_Total_FObj(k,:)));
    end
    FObj_Indiv2(j) = (Index1 + Index2 + Index3 + 1);
  end
  
  // We compute the rank for Pop
  for j=1:size(Aux_Total_FObj,1)
    Index1 = 0; Index2 = 0; Index3 = 0;
    for k=1:size(Aux_FObj_Indiv1,1)
      Index1 = Index1 + double(and(Aux_Total_FObj(j,:)<=Aux_FObj_Indiv1(k,:)) & or(Aux_Total_FObj(j,:)<Aux_FObj_Indiv1(k,:)));
      Index2 = Index2 + double(and(Aux_Total_FObj(j,:)<=Aux_FObj_Indiv2(k,:)) & or(Aux_Total_FObj(j,:)<Aux_FObj_Indiv2(k,:)));
    end
    for k=1:size(Aux_Total_FObj,1)
      Index3 = Index3 + double(and(Aux_Total_FObj(j,:)<=Aux_Total_FObj(k,:)) & or(Aux_Total_FObj(j,:)<Aux_Total_FObj(k,:)));
    end
    Total_FObj(j) = (Index1 + Index2 + Index3 + 1);
  end
  //
  // Recombination
  //
  select selection_strategy
  case 'elitist' then
    // Fusion of the various populations
    Total_Pop  = lstcat(Total_Pop, Indiv1, Indiv2);
    Total_FObj = [Total_FObj' FObj_Indiv1' FObj_Indiv2']';  
    Aux_Total_FObj = [Aux_Total_FObj' Aux_FObj_Indiv1' Aux_FObj_Indiv2']';  
    // Normalization of the efficiency
    FObj_Pop_Max = max(Total_FObj);
    FObj_Pop_Min = min(Total_FObj);
    
    Efficiency = (1 - pressure) * (Total_FObj - FObj_Pop_Min)/max([FObj_Pop_Max - FObj_Pop_Min, %eps]) + pressure;
    [Efficiency, Index_Sort] = sort(Efficiency);
    Efficiency = Efficiency(1:pop_size);
    // Extraction and selection of the phenotype
    // Total_FObj
    Total_FObj     = Total_FObj(Index_Sort);
    Total_FObj     = Total_FObj(1:pop_size);
    // Total_FObj
    Aux_Total_FObj = Aux_Total_FObj(Index_Sort,:);
    Aux_Total_FObj = Aux_Total_FObj(1:pop_size,:);
    // Extraction and selection of the genotype
    Total_Pop = list(Total_Pop(Index_Sort));
    Total_Pop = list(Total_Pop(1:pop_size));
    if (Log) then
      printf(' - min / max value found = %f / %f\n', FObj_Pop_Min, FObj_Pop_Max);
    end
  else
    error('optim_ga: wrong selection strategy');
  end
end

pop_opt      = codage_func(Total_Pop,'decode');
fobj_pop_opt = Aux_Total_FObj;
endfunction

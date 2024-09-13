function [pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_nsga2(ga_f, pop_size, nb_generation, p_mut, p_cross, init_func, param, Log, crossover_func, mutation_func, codage_func, nb_couples)

[nargout, nargin] = argn();

if ~isdef('ga_f','local') then
  error('optim_nsga2: ga_f is mandatory');
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
if ~isdef('nb_couples','local') then
  nb_couples = 100;
end

// Initialization of the population
if (Log) then
  printf('optim_nsga2: Initialization of the population\n');
end

Total_Pop = init_func(pop_size);
// Code the individuals
Total_Pop = codage_func(Total_Pop,'code');

for i=1:length(Total_Pop)
  Total_FObj(i,:) = ga_f(Total_Pop(i));
end
Aux_Total_FObj = Total_FObj;

// Compute the domination rank
Index = 1:size(Aux_Total_FObj,1);
Rank  = zeros(size(Aux_Total_FObj,1),1);
Count = 1;
while size(Aux_Total_FObj,1)>1
  [tmp1,tmp2,Index_List]  = pareto_filter(Aux_Total_FObj);
  Rank(Index(Index_List)) = Count;
  Count = Count + 1;
  Aux_Total_FObj(Index_List,:) = [];
  Index(Index_List) = [];
end

// Compute the crowding distance 
Aux_Total_FObj = Total_FObj;

Index    = 1:size(Aux_Total_FObj,1);
Crowdist = zeros(size(Aux_Total_FObj,1),1);
for i=1:size(Total_FObj,2)
  [tmp, Index_List] = sort(Aux_Total_FObj(:,i));
  Aux_Total_FObj    = Aux_Total_FObj(Index_List,:);
  Index             = Index(Index_List);
  Crowdist(Index_List(1)) = %inf;
  Crowdist(Index_List($)) = %inf;
  _Max = max(Aux_Total_FObj(:,i));
  _Min = min(Aux_Total_FObj(:,i));
  for j=2:size(Aux_Total_FObj,1)-1
    Crowdist(Index(j)) = Crowdist(Index(j)) - (Aux_Total_FObj(j+1,i) - Aux_Total_FObj(j-1,i)) / (_Max - _Min);
  end
end

if (nargout==4) then
  pop_init = Total_Pop;
  fobj_pop_init = Total_FObj;
end

// The genetic algorithm
for It=1:nb_generation
  if (Log) then
    printf('optim_nsga2: iteration %d / %d\n', It, nb_generation);
  end
  //
  // Selection
  //
  Indiv1 = list();
  Indiv2 = list();
  for j=1:nb_couples
    // Selection of 2 individuals via binary tournament selection to fill Indiv1
    Index1 = ceil((size(Total_FObj,1) - 1)*rand(1,1)+1);
    Index2 = ceil((size(Total_FObj,1) - 1)*rand(1,1)+1);
    if (Rank(Index1)<Rank(Index2)) | ((Rank(Index1)==Rank(Index2)) & (Crowdist(Index1)>Crowdist(Index2))) then
      Indiv1(j) = Total_Pop(Index1);
    else
      Indiv1(j) = Total_Pop(Index2);
    end
    // Selection of 2 individuals via binary tournament selection to fill Indiv2
    Index1 = ceil((size(Total_FObj,1) - 1)*rand(1,1)+1);
    Index2 = ceil((size(Total_FObj,1) - 1)*rand(1,1)+1);
    if (Rank(Index1)<Rank(Index2)) | ((Rank(Index1)==Rank(Index2)) & (Crowdist(Index1)>Crowdist(Index2))) then
      Indiv2(j) = Total_Pop(Index1);
    else
      Indiv2(j) = Total_Pop(Index2);
    end
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
    FObj_Indiv1(j,:) = ga_f(Indiv1(j)); 
    FObj_Indiv2(j,:) = ga_f(Indiv2(j));
  end

  // We merge all the individuals in one list ...  
  All_Pop  = lstcat(Total_Pop, Indiv1, Indiv2);
  All_FObj = [Total_FObj' FObj_Indiv1' FObj_Indiv2']';  

  Aux_All_FObj = All_FObj;

  // Compute the domination rank on all the population
  Index = 1:size(Aux_All_FObj,1);
  Rank  = zeros(size(Aux_All_FObj,1),1);
  Count = 1;
  while size(Aux_All_FObj,1)>1
    [tmp1,tmp2,Index_List]  = pareto_filter(Aux_All_FObj);
    Rank(Index(Index_List)) = Count;
    Count = Count + 1;
    Aux_All_FObj(Index_List,:) = [];
    Index(Index_List)          = [];
  end

  // Compute the crowding distance
  Aux_All_FObj = All_FObj;
  
  Index    = 1:size(Aux_All_FObj,1);
  Crowdist = zeros(size(Aux_All_FObj,1),1);
  for k=1:size(Aux_All_FObj,2)
    [tmp, Index_List] = sort(Aux_All_FObj(:,k));
    Aux_All_FObj = Aux_All_FObj(Index_List,:);
    Index = Index(Index_List);
    Crowdist(Index_List(1)) = %inf;
    Crowdist(Index_List($)) = %inf;
    _Max = max(Aux_All_FObj(:,k));
    _Min = min(Aux_All_FObj(:,k));
    for j=2:size(Aux_All_FObj,1)-1
      Crowdist(Index(j)) = Crowdist(Index(j)) - (Aux_All_FObj(j+1,k) - Aux_All_FObj(j-1,k)) / (_Max - _Min);
    end
  end
  //
  // Recombination
  //
  // We rank all the individual wrt to the partial order
  for k=1:size(All_FObj,1)-1
    for j=k+1:size(All_FObj,1)
      if (Rank(j)<Rank(k)) | ((Rank(j)==Rank(k)) & (Crowdist(j)>Crowdist(k))) then
        tmp           = Rank(k);
        Rank(k)       = Rank(j);
        Rank(j)       = tmp;
        tmp           = Crowdist(k);
        Crowdist(k)   = Crowdist(j);
        Crowdist(j)   = tmp;
        tmp           = All_Pop(k);
        All_Pop(k)    = All_Pop(j);
        All_Pop(j)    = tmp;
        tmp           = All_FObj(k,:);
        All_FObj(k,:) = All_FObj(j,:);
        All_FObj(j,:) = tmp;
      end
    end
  end
  // Extraction and selection of the phenotype
  Total_FObj = All_FObj(1:pop_size,:);
  // Extraction and selection of the genotype
  Total_Pop = list(All_Pop(1:pop_size));
  // Extraction of the ranks and Crow distance
  Rank     = Rank(1:pop_size);
  Crowdist = Crowdist(1:pop_size);
end

pop_opt      = codage_func(Total_Pop,'decode');
fobj_pop_opt = Total_FObj;
endfunction

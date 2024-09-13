// example of use of the genetic algorithm

func = 'rosenbrock';

PopSize     = 100;
Proba_cross = 0.7;
Proba_mut   = 0.1;
NbGen       = 5;
NbCouples   = 110;
Strategy    = 'elitist';
Log         = %T;
nb_disp     = 10; // Nb point to display from the optimal population
pressure    = 0.05;

ga_params = init_param();
ga_params = add_param(ga_params,'minbound',eval('get_min_bound_'+func+'()'));
ga_params = add_param(ga_params,'maxbound',eval('get_max_bound_'+func+'()'));
ga_params = add_param(ga_params,'dimension',2);
ga_params = add_param(ga_params,'beta',0);
ga_params = add_param(ga_params,'delta',0.1);

//////////////////////////////////////////

drawlater;
Min = get_param(ga_params,'minbound');
Max = get_param(ga_params,'maxbound');
x = Min(1):(Max(1)-Min(1))/20:Max(1); y = Min(2):(Max(2)-Min(2))/20:Max(2);
[X,Y]=meshgrid(x,y);
for i=1:size(X,1)
  for j=1:size(X,2)
    Z(i,j) = eval(func+'([X(i,j) Y(i,j)])');
  end
end
xset('fpf',' ');
contour(x,y,Z', 10);
_axes = get("current_axes");
_axes.data_bounds = [Min(1) Max(1) Min(2) Max(2)];
deff('y=f(x)','y = '+func+'(x)');
x0  = (Max - Min) .* rand(size(Min,1),size(Min,2)) + Min;

drawnow;

///////////////////////
// Genetic Algorithm //
///////////////////////

[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init] = optim_ga(f, PopSize, NbGen, Proba_mut, Proba_cross, init_func_default, ga_params, Log, ...
                                                            crossover_func_default, mutation_func_default, codage_identity, Strategy, NbCouples, pressure);

if (size(pop_opt(1),2)==2) then
  drawlater;
  printf('plotting init population ...\n');
  for i=1:length(pop_init)
    plot(pop_init(i)(1),pop_init(i)(2),'r.');
  end
  printf('plotting result population ...\n');
  for i=1:length(pop_opt)
    plot(pop_opt(i)(1),pop_opt(i)(2),'g.');
  end
  drawnow;
end

printf('Genetic Algorithm: %d points from pop_opt\n', nb_disp); 
disp(pop_opt(1:nb_disp))
printf('Genetic Algorithm: %d points from fobj_pop_opt\n', nb_disp); 
disp(fobj_pop_opt(1:nb_disp))

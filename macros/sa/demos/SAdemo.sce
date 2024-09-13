// example of use of the simulated annealing method

func = 'rbf'; // ell = quadratic
              // rbf = Rosenbrock
              // prbl = parabol
              // cuter = CUTEr function
CuterFunction = 'PALMER1C'; // Name of the CUTEr function used (a list of function can be obtained via
                            // the command cuter_get_list_of_pbs_unconstrained)

Proba_start = 0.8;
It_intern = 1000;
It_extern = 30;
It_Pre    = 100;

//////////////////////////////////////////

printf('\n*************** start example %s ***************\n\n', func);
drawlater;
select func
  case 'rbf' then
    x=-1.5:0.01:2; y=-0.5:0.01:2;
    [X,Y]=meshgrid(x,y);
    Z=rbf(X,Y);
    LZ=log(1.0e-16+Z);
    contour(x,y,LZ', log([500 200 50 10 4 1 0.1]));
    _axes = get("current_axes");
    _axes.data_bounds = [-1.5 2 -0.5 2];
    f   = rbf;
    x0  = [-1.2;1];
  case 'ell' then
    x=-1:0.005:1; y=-1:0.005:1;
    [X,Y]=meshgrid(x,y);
    Z=ell(X,Y);
    contour(x,y,Z', [4 3 2 1 0.5 0.25 0.1 0.01 0.001 0.0001]);
    _axes = get("current_axes");
    _axes.data_bounds = [-1 1 -1 1];
    f   = ell;
    x0  = [0.5;0.25];
  case 'prbl' then
    x = -2:0.01:2; y = -2:0.01:2;
    [X,Y]=meshgrid(x,y);
    Z  = parabol(X,Y);
    contour(x,y,Z', 10);
    _axes = get("current_axes");
    _axes.data_bounds = [-2 2 -2 2];
    f   = parabol;
    x0  = [1;1];
  case 'cuter' then
    f   = cuter_fobj;
    [x0, lb, ub] = cuter_init_unconstrained(CuterFunction);
  else
    printf('unknown function %s', func);
end
drawnow;

/////////////////////////
// Simulated Annealing //
/////////////////////////

printf('SA: geometrical decrease temperature law\n');

T0 = compute_initial_temp(x0, f, Proba_start, It_Pre, neigh_func_default);
printf('Initial temperatore T0 = %f\n', T0);

[x_opt, f_opt, sa_mean_list, sa_var_list, temp_list] = optim_sa(x0, f, It_extern, It_intern, T0, Log = %T);

printf('optimal solution:\n'); disp(x_opt);
printf('value of the objective function = %f\n', f_opt);

scf();
drawlater;
subplot(2,1,1);
xtitle('Geometrical annealing','Iteration','Mean / Variance');
t = 1:length(sa_mean_list);
plot(t,sa_mean_list,'r',t,sa_var_list,'g');
legend(['Mean','Variance']);
subplot(2,1,2);
xtitle('Temperature evolution','Iteration','Temperature');
for i=1:length(t)-1
  plot([t(i) t(i+1)], [temp_list(i) temp_list(i)],'k-');
end
drawnow;

/////////
// FSA //
/////////

printf('SA: the FSA algorithm\n');

T0 = compute_initial_temp(x0, f, Proba_start, It_Pre, neigh_func_default);
printf('Initial temperatore T0 = %f\n', T0);

[x_opt, f_opt, sa_mean_list, sa_var_list, temp_list] = optim_sa(x0, f, It_extern, It_intern, T0, Log = %T, ...
                                                                temp_law_fsa,   [], ...
                                                                neigh_func_fsa, []);

printf('optimal solution:\n'); disp(x_opt);
printf('value of the objective function = %f\n', f_opt);

scf();
drawlater;
subplot(2,1,1);
xtitle('FSA','Iteration','Mean / Variance');
t = 1:length(sa_mean_list);
plot(t,sa_mean_list,'r',t,sa_var_list,'g');
legend(['Mean','Variance']);
subplot(2,1,2);
xtitle('Temperature evolution','Iteration','Temperature');
for i=1:length(t)-1
  plot([t(i) t(i+1)], [temp_list(i) temp_list(i)],'k-');
end
drawnow;

//////////
// VFSA //
//////////

printf('SA: the VFSA algorithm\n');

T0 = compute_initial_temp(x0, f, Proba_start, It_Pre, neigh_func_default);
printf('Initial temperatore T0 = %f\n', T0);

[x_opt, f_opt, sa_mean_list, sa_var_list, temp_list] = optim_sa(x0, f, It_extern, It_intern, T0, Log = %T, ...
                                                                temp_law_vfsa,[], ...
                                                                neigh_func_vfsa, []);

printf('optimal solution:\n'); disp(x_opt);
printf('value of the objective function = %f\n', f_opt);

scf();
drawlater;
subplot(2,1,1);
xtitle('VFSA','Iteration','Mean / Variance');
t = 1:length(sa_mean_list);
plot(t,sa_mean_list,'r',t,sa_var_list,'g');
legend(['Mean','Variance']);
subplot(2,1,2);
xtitle('Temperature evolution','Iteration','Temperature');
for i=1:length(t)-1
  plot([t(i) t(i+1)], [temp_list(i) temp_list(i)],'k-');
end
drawnow;

/////////
// CSA //
/////////

printf('SA: the CSA algorithm\n');

T0 = compute_initial_temp(x0, f, Proba_start, It_Pre, neigh_func_default);
printf('Initial temperatore T0 = %f\n', T0);

[x_opt, f_opt, sa_mean_list, sa_var_list, temp_list] = optim_sa(x0, f, It_extern, It_intern, T0, Log = %T, ...
                                                                temp_law_csa, [], ...
                                                                neigh_func_csa, []);

printf('optimal solution:\n'); disp(x_opt);
printf('value of the objective function = %f\n', f_opt);

scf();
drawlater;
subplot(2,1,1);
xtitle('Classical simulated annealing','Iteration','Mean / Variance');
t = 1:length(sa_mean_list);
plot(t,sa_mean_list,'r',t,sa_var_list,'g');
legend(['Mean','Variance']);
subplot(2,1,2);
xtitle('Temperature evolution','Iteration','Temperature');
for i=1:length(t)-1
  plot([t(i) t(i+1)], [temp_list(i) temp_list(i)],'k-');
end
drawnow;

///////////
// Huang //
///////////

printf('SA: the Huang annealing\n');

T0 = compute_initial_temp(x0, f, Proba_start, It_Pre, neigh_func_default);
printf('Initial temperatore T0 = %f\n', T0);

[x_opt, f_opt, sa_mean_list, sa_var_list, temp_list] = optim_sa(x0, f, It_extern, It_intern, T0, Log = %T, ...
                                                                temp_law_huang, [], ... 
                                                                neigh_func_default, []);

printf('optimal solution:\n'); disp(x_opt);
printf('value of the objective function = %f\n', f_opt);

scf();
drawlater;
subplot(2,1,1);
xtitle('Huang annealing','Iteration','Mean / Variance');
t = 1:length(sa_mean_list);
plot(t,sa_mean_list,'r',t,sa_var_list,'g');
legend(['Mean','Variance']);
subplot(2,1,2);
xtitle('Temperature evolution','Iteration','Temperature');
for i=1:length(t)-1
  plot([t(i) t(i+1)], [temp_list(i) temp_list(i)],'k-');
end
drawnow;


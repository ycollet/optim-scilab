// DSSAdemo.sce

func     = 'rbf'; // ell   = quadratic
                  // rbf   = Rosenbrock
                  // prbl  = parabol
                  // cuter = CUTEr function
ItMX          = 150;  // maximum number of descent steps

CuterFunction = 'PALMER1C'; // Name of the CUTEr function used (a list of function can be obtained via
                            // the command cuter_get_list_of_pbs_unconstrained)

lines(0);

rand('seed',12345);

//////////////////////////////////////////

printf('\n*************** start example %s ***************\n\n', func);

clf;
drawlater;

select func
  case 'rbf' then
    x=-1.5:0.01:2; y=-0.5:0.01:2;
    [X,Y]=meshgrid(x,y);
    Z=rbf(X,Y);
    Z=log(1.0e-16+Z);
    contour(x,y,Z', log([500 200 50 10 4 1 0.1]));
    _axes = get("current_axes");
    _axes.data_bounds = [-1.5 2 -0.5 2];
    f   = rbf;
    x0  = [-1.1;1];
    // For plot_d
    lbounds = [-1.5;-0.5];
    ubounds = [2;2];
  case 'ell' then
    x=-1:0.005:1; y=-1:0.005:1;
    [X,Y]=meshgrid(x,y);
    Z=ell(X,Y);
    contour(x,y,Z', [4 3 2 1 0.5 0.25 0.1 0.01 0.001 0.0001]);
    _axes = get("current_axes");
    _axes.data_bounds = [-1 1 -1 1];
    f   = ell;
    x0  = [0.5;0.25];
    // For plot_d
    lbounds = [-1;-1];
    ubounds = [1;1];
  case 'prbl' then
    x = -2:0.01:2; y = -2:0.01:2;
    [X,Y]=meshgrid(x,y);
    Z  = parabol(X,Y);
    contour(x,y,Z', 10);
    _axes = get("current_axes");
    _axes.data_bounds = [-2 2 -2 2];
    f   = parabol;
    x0  = [1;1];
    // For plot_d
    lbounds = [-2;-2];
    ubounds = [2;2];
  case 'cuter' then
    f   = cuter_fobj;
    [x0, lb, ub] = cuter_init_unconstrained(CuterFunction);
    // For plot_d
    lbounds = lb;
    ubounds = ub;
  else
    printf('unknown function %s', func);
end

xgrid(2);

////////////////////// Nelder and Mead optimization

[x_opt, x_history] = optim_DSSA(f, ubounds, lbounds, ItMX);

printf('Number of function evaluation: %d\n',length(x_history));
printf('dimension 1 of x_history = %d\n',length(x_history(1)));

printf('xopt = ');
disp(x_opt)
printf('f_opt = %f\n', f(x_opt));

drawlater;
for i=4:length(x_history)
  plot(x_history(i)(1)(1), x_history(i)(1)(2), 'ko');
  plot(x_history(i)(2)(1), x_history(i)(2)(2), 'ko');
  plot(x_history(i)(3)(1), x_history(i)(3)(2), 'ko');
  plot([x_history(i)(1)(1) x_history(i)(2)(1)], [x_history(i)(1)(2) x_history(i)(2)(2)], 'k-');
  plot([x_history(i)(2)(1) x_history(i)(3)(1)], [x_history(i)(2)(2) x_history(i)(3)(2)], 'k-');
  plot([x_history(i)(3)(1) x_history(i)(1)(1)], [x_history(i)(3)(2) x_history(i)(1)(2)], 'k-');
end
xtitle('Evolution of the simplex','x1','x2');
drawnow;


func         = 'rbf';  // ell = quadratic, rbf = Rosenbrock, prbl = parabol
cgm          = 'rosen0'; // rosen0 = original rosenbrock method
                         // rosen1 = rosenbrock method with swann modification
                         // rosen2 = rosenbrock method with palmer modification
ItMX         = 3000;     // maximum number of descent steps
TOL          = 1e-4;
XTOL         = 1.0e-4; // accuracy for convergence test (minimum)
StepTOL      = 1.0e-4; // minimal size of the step before restarting
Log          = %T;      // print info from line step, 1 = on / 0 = off
_alpha       = 2;    // Increase step size coefficient
_beta        = 0.5;  // Decrease step size coefficient
h            = 0.01; // step size


///////////////////////////////////////

printf('\n*************** start example %s ***************\n\n', func);

clf;
drawlater;

select func
  case 'rbf' then
    x=-1.5:0.01:2; y=-0.5:0.01:2;
    [X,Y]=meshgrid(x,y);
    Z=rbf(X,Y);
    LZ=log(1.0e-16+Z);
    contour(x,y,LZ', log([500 200 50 10 4 1 0.1]))
    _axes = get("current_axes");
    _axes.data_bounds = [-1.5 2 -0.5 2];
    f   = rbf;
    fs  = grad_rbf;
    fss = hess_rbf;
    x0  = [-1.2;1];
   case 'ell' then
     x=-1:0.005:1; y=-1:0.005:1;
    [X,Y]=meshgrid(x,y);
    Z=ell(X,Y);
    contour(x,y,Z', [4 3 2 1 0.5 0.25 0.1 0.01 0.001 0.0001]);
    _axes = get("current_axes");
    _axes.data_bounds = [-1 1 -1 1];
    f   = ell;
    fs  = grad_ell;
    fss = hess_ell;
    x0  = [0.5;0.25];
  case 'prbl' then
    x = -2:0.01:2; y = -2:0.01:2;
    [X,Y]=meshgrid(x,y);
    Z  = parabol(X,Y);
    contour(x,y,Z', 10);
    _axes = get("current_axes");
    _axes.data_bounds = [-2 2 -2 2];
    f   = parabol;
    fs  = grad_parabol;
    fss = hess_parabol;
    x0  = [1;1];
  else
    printf('unknown function %s', func);
end

xgrid(2);

////////////////// Rosenbrock method

[x_opt, x_history] = optim_rosen(x0, f, h, Log, ItMX, TOL, StepTOL, XTOL, cgm);

printf('Number of function evaluation: %d\n',length(x_history));
printf('Final solution = %f\n', f(x_opt)); disp(x_opt');

plot(x_history(1)(1), x_history(1)(2), 'gx');
for i=2:length(x_history)-1
  plot([x_history(i)(1) x_history(i-1)(1)], [x_history(i)(2) x_history(i-1)(2)], 'k-');
  plot(x_history(i)(1), x_history(i)(2), 'ro');
end
n = length(x_history);
plot([x_history(n)(1) x_history(n-1)(1)], [x_history(n)(2) x_history(n-1)(2)], 'k-');
plot(x_history(n)(1), x_history(n)(2), 'gx');

xtitle('La méthode de Rosenbrock','x1','x2');
drawnow;

printf('powelldemo: minimum not found, last point:\nx=[%.4e %.4e], ||f''(x)||=%4e\n', x_opt(1), x_opt(2), fs(x_opt));

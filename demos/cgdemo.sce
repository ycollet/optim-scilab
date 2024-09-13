func         = 'rbf'; // ell   = quadratic
                        // rbf   = Rosenbrock
                        // prbl  = parabol
                        // cuter = CUTEr function
cgm          = 'pr';   // fr = fletcher-reeves
                       // pr = polak-ribiere
                       // hs = hestenes-stiefel
                       // pw = powel
                       // dy = Dai-Yuan
ItMX         = 100;    // maximum number of descent steps
ls_ItMX      = 10;    // maximum number of line search steps
Restart      = 10;     // restart cg every Restart iterations
TOL          = 1.0e-6; // accuracy for convergence test (minimum)
StepTOL      = 1.0e-6; // accuracy for convergence test - size of the step 
XTOL         = 1.0e-6; // accuracy for convergence test - improvement on x between two iterations
SubTOL       = 1.0e-6; // accuracy for convergence test (line search)
Log          = %T;      // print info from line step, 1 = on / 0 = off
h            = 0.01;  // initial step length in line search
PowRestart   = %F;     // Powell restart (we force a powell restart at the first iteration)
nu           = 0.1;    // level of similitude before powell restart
StepIncrease = 0.01;
CuterFunction = 'PALMER3C'; // Name of the CUTEr function used (a list of function can be obtained via
                            // the command cuter_get_list_of_pbs_unconstrained)
Plot         = %T;
Spectral     = %F;
Restart      = 7;

//lsm  = ls_secant;
//lsm  = ls_dicho;
//lsm  = ls_goldsect;
lsm  = ls_backtrack;
//lsm  = ls_polynom;
///////////////////////////////////////

printf('\n*************** start example %s ***************\n\n', func);

clf;

lines(0);

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
    fs  = grad_ell;
    fss = hess_ell;
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
    fs  = grad_parabol;
    fss = hess_parabol;
    x0  = [1;1];
    // For plot_d
    lbounds = [-2;-2];
    ubounds = [2;2];
  case 'cuter' then
    f   = cuter_fobj;
    fs  = cuter_grad;
    fss = cuter_hess;
    [x0, lb, ub] = cuter_init_unconstrained(CuterFunction);
    // For plot_d
    lbounds = lb;
    ubounds = ub;
  else
    printf('unknown function %s', func);
end

xgrid(2);

////////////////// conjugate gradients

printf('cgdemo: start point:\nx=[%.4e %.4e], ||f''(x)||=%4e\n', x0(1), x0(2), fs(x0));

[x_opt, x_history] = optim_cg(x0, f, fs, h, Log, ItMX, ls_ItMX, lsm, cgm, TOL, StepTOL, XTOL, SubTOL, Restart, Spectral);

printf('cgdemo: minimum not found, last point:\nx=[%.4e %.4e], ||f''(x)||=%4e\n', x_opt(1), x_opt(2), fs(x_opt));

if (Plot) then
  drawlater;
  for i=1:length(x_history)
    plot([x_history(i)(1)(1) x_history(i)(1)(1)], [x_history(i)(1)(2) x_history(i)(1)(2)], 'ro');
    if (length(x_history(i))>1) then
      for j=2:length(x_history(i))
        plot([x_history(i)(j)(1) x_history(i)(j)(1)], [x_history(i)(j)(2) x_history(i)(j)(2)], 'gx');
        plot([x_history(i)(j-1)(1) x_history(i)(j)(1)], [x_history(i)(j-1)(2) x_history(i)(j)(2)], 'k-');
      end
    end
  end
  plot([x_opt(1) x_opt(1)], [x_opt(2) x_opt(2)], 'ro');
  xtitle('Conjugate gradient - '+cgm,'x1','x2');
  drawnow;
end

NbIter = 0;
for i=1:length(x_history)
  NbIter = NbIter + length(x_history(i));
end
printf('Number of call to the objective function = %d\n',NbIter);

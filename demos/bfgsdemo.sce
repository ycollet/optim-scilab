//funcname = 'rosenbrock';
funcname = 'mccormic';
//funcname = 'sixhumpcamelback';
//funcname = 'branin2'; // BFGS doesn't work
//funcname = 'schubert';// BFGS doesn't work
//funcname = 'hansen';// BFGS doesn't work
//funcname = 'paviani';
//funcname = 'booth';
//funcname = 'matyas';
//funcname = 'rastrigin';
//funcname = 'griewank2'; // BFGS does't work
//funcname = 'exp2'; // Infinite loop
//funcname = 'treccani';
//funcname = 'branin';
//funcname = 'colville';
//funcname = 'chichinadze';
//funcname = 'hartmann34';
//funcname = 'hartmann64';
//funcname = 'price';
//funcname = 'goldsteinprice';
//funcname = 'dixonprice';
//funcname = 'hump';
//funcname = 'dejongf2';
//funcname = 'dejongf5';// Division by zero - ls_secant
//funcname = 'dejongf7';
//funcname = 'schafferf6'
//funcname = 'schafferf7';
//funcname = 'stuckman'; // Division by zero - ls_secant
//funcname = 'easom'; // Division by zero - ls_secant
//funcname = 'bohachevsky1';
//funcname = 'bohachevsky2';
//funcname = 'bohachevsky3';
//funcname = 'beale';
//funcname = 'levy13'; // BFGS doesn't work
//funcname = 'levy8'; // BFGS doesn't work
//funcname = 'levy5';
//funcname = 'levy2'; // BFGS doesn't work
//funcname = 'holtzmann';
//funcname = 'gen_rosen';
//funcname = 'shekel';
//funcname = 'griewank';
//funcname = 'sphere';
//funcname = 'weierstrass'; // BFGS doesn't work
//funcname = 'ackley';
//funcname = 'ellipsoid';
//funcname = 'rotell';
//funcname = 'abspow';
//funcname = 'michalewicz';
//funcname = 'powell'; // Division by zero - ls_secant
//funcname = 'power'; // Long - probably doesn't work
//funcname = 'gen_rastrigin';
//funcname = 'schwefel'; // Division by zero - ls_secant - sometime it works ||df||=0
//funcname = 'trid';
//funcname = 'zhakarov';
//funcname = 'freudroth';
//funcname = 'himmelblau';
//funcname = 'jensamp';
//funcname = 'zhufu'; // BFGS doen't work
//funcname = 'cola';
//funcname = 'leon';
//funcname = 'giunta';
//funcname = 'bukin2';
//funcname = 'bukin4'; // BFGS doesn't work - works with ls_dicho
//funcname = 'bukin6'; // BFGS doesn't work - initial direction not a descent direction
//funcname = 'stybtang';
//funcname = 'zettl';
//funcname = 'threehumpcamelback';

Max = eval('get_max_bound_' + funcname + '(2)');
Min = eval('get_min_bound_' + funcname + '(2)');

deff('y=fobjs(x)','y = ' + funcname + '(x);');
deff('y=dfobjs(x)','y = derivative(' + funcname + ',x)'';');

x0 = (Max-Min).*rand(size(Min,1),size(Min,2)) + Min;

ItMX     = 500;         // maximum number of descent steps
ls_ItMX  = 100;        // maximum number of line search iteration
Restart  = 10;         // restart cg every Restart iterations
TOL      = 1.0e-6;     // accuracy for convergence test (minimum)
StepTOL  = 1.0e-12;     // accuracy for convergence test - size of the step 
XTOL     = 1.0e-12;     // accuracy for convergence test - improvement on x between two iterations
SubTOL   = 1.0e-4;     // accuracy for convergence test (line search)
Log      = %F;          // print info from line step, 1 = on / 0 = off
h        = 4.0; //0.01;       // initial step length in line search

Plot     = %T;

// lsm  = ls_secant;
lsm  = ls_dicho;
// lsm  = ls_goldsect;
// lsm  = ls_backtrack;
// lsm  = ls_polynom;

////////////// Broyden Fletcher Goldfarb Shanno quasi-Newton

[x_opt, x_history] = optim_bfgs(x0, fobjs, dfobjs, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL);

printf('bfgsdemo: last point: ||f''(x)||=%4e\n', norm(dfobjs(x_opt))); 
printf('bfgsdemo: x-opt = '); disp(x_opt');

if (Plot) then
  drawlater;
  X = Min(1):(Max(1)-Min(1))/10:Max(1);
  Y = Min(2):(Max(2)-Min(2))/10:Max(2);
  for i=1:length(X)
    for j=1:length(Y)
      Z(i,j) = eval(funcname+'([X(i),Y(j)])');
    end
  end
  contour(X,Y,Z,10);
  
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
  xtitle('BFGS Newton','x1','x2');
  drawnow;
end

NbIter = 0;
for i=1:length(x_history)
  NbIter = NbIter + length(x_history(i));
end
printf('Number of call to the objective function = %d\n',NbIter);

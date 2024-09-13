functoplot='rosenbrock';
//functoplot='mccormic';
//functoplot='sixhumpcamelback';
//functoplot='branin';
//functoplot='schubert';
//functoplot='hansen';
//functoplot='paviani';
//functoplot='booth';
//functoplot='matyas';
//functoplot='rastrigin';
//functoplot='griewank2';
//functoplot='exp2';
//functoplot='treccani';
//functoplot='branin2';
//functoplot='chichinadze';
//functoplot='price';
//functoplot='goldsteinprice';
//functoplot='dejongf2';
//functoplot='dejongf5';
//functoplot='dejongf7';
//functoplot='schafferf6';
//functoplot='schafferf7';
//functoplot='stuckman';
//functoplot='easom';
//functoplot='bohachevsky1';
//functoplot='bohachevsky2';
//functoplot='beale';
//functoplot='holtzmann';
//functoplot='griewank';
//functoplot='sphere';
//functoplot='weierstrass';
//functoplot='ackley';
//functoplot='ellipsoid';
//functoplot='rotell';
//functoplot='abspow';
//functoplot='michalewicz';
//functoplot='levy5';
//functoplot='levy13';

deff('y=fobj(x)','y='+functoplot+'(x)');
deff('y=dfobj(x)','[tmp,y]=CDcost(x,3,'+functoplot+'); y=y;'); // compute only df

// example of use of the simulated tempering method

ItMX     = 50;         // maximum number of descent steps
ls_ItMX  = 100;        // maximum number of line search iteration
Restart  = 10;         // restart cg every Restart iterations
TOL      = 1.0e-6;     // accuracy for convergence test (minimum)
StepTOL  = 1.0e-12;    // accuracy for convergence test - size of the step 
XTOL     = 1.0e-12;    // accuracy for convergence test - improvement on x between two iterations
SubTOL   = 1.0e-4;     // accuracy for convergence test (line search)
Log      = %F;         // print info from line step, 1 = on / 0 = off
h        = 4.0; //0.01;       // initial step length in line search

// lsm  = ls_secant;
lsm  = ls_dicho;
// lsm  = ls_goldsect;
// lsm  = ls_backtrack;
// lsm  = ls_polynom;

//////////////////////////////////////////

x_min = eval('get_min_bound_'+functoplot+'()');
x_max = eval('get_max_bound_'+functoplot+'()');

x0 = (x_max - x_min) .* rand(size(x_min,1),size(x_min,2)) + x_min;
x0 = [-1.2; 1];
x = x_min(1):(x_max(1) - x_min(1))/ 10:x_max(1);
y = x_min(2):(x_max(2) - x_min(2))/ 10:x_max(2);

[X,Y]=meshgrid(x,y);
for i=1:size(X,1)
  for j=1:size(X,2)
    Z(i,j) = eval(functoplot+'([X(i,j) Y(i,j)])');
  end
end

//////////
// BFGS //
//////////

printf('BFGS: starting optimization\n');

[x_opt, x_history] = optim_bfgs(x0, fobj, dfobj, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL);

printf('Solution found:'); disp(x_opt');
printf('value of the objective function = %f\n', fobj(x_opt));
printf('Optimal solution:'); disp(eval('get_opti_'+functoplot+'()')');

drawlater;

xset('fpf',' ');
contour(x,y,Z,10);
_axes = get("current_axes");
_axes.data_bounds = [x_min(1) x_max(1) x_min(2) x_max(2)];

for i=1:length(x_history)
  for j=1:length(x_history(i))-1
    plot([x_history(i)(j)(1) x_history(i)(j+1)(1)],[x_history(i)(j)(2) x_history(i)(j+1)(2)],'k-');
  end
  plot(x_history(i)($)(1),x_history(i)($)(2),'ro');
end
drawnow;


//functoplot='rosenbrock';
//functoplot='mccormic';
//functoplot='sixhumpcamelback';
//functoplot='branin';
functoplot='schubert';
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
deff('y=dfobj(x)','y=derivative(fobj,x)'';');

ItMX     = 50;         // maximum number of descent steps
TOL      = 1.0e-6;     // accuracy for convergence test - derivatives 
XTOL     = 1.0e-6;     // accuracy for convergence test - improvement on x between two iterations
Log      = %T;          // print info from line step, 1 = on / 0 = off
DXMax    = 1.1;

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

//////////////////////
// Barzilai-Borwein //
//////////////////////

[x_opt, x_history] = optim_bb(x0, f, fs, Log, ItMX, DXMax, TOL, XTOL);

printf('bbdemo: starting point - fobj = %f - x =', fobj(x0)); disp(x0');

printf('bbdemo: last point - fobj = %f - x = :', fobj(x_opt)); disp(x_opt');

x_opti = eval('get_opti_'+functoplot+'()');
if ~isnan(x_opti(1)) then
  f_opti = fobj(x_opti(:,1));
else
  f_opti = %nan;
end
printf('bbdemo: optimum for the problem %s - fobj = %f - x = :',functoplot, f_opti); disp(x_opti');

drawlater;
xset('fpf',' ');
contour(x,y,Z,10);
_axes = get("current_axes");
_axes.data_bounds = [x_min(1) x_max(1) x_min(2) x_max(2)];

for i=1:length(x_history)-1
  plot([x_history(i)(1) x_history(i+1)(1)], [x_history(i)(2) x_history(i+1)(2)], 'k-');
  plot(x_history(i)(1), x_history(i)(2), 'ro');
end
plot([x_opt(1) x_opt(1)], [x_opt(2) x_opt(2)], 'ro');
xtitle('Barzilai Borwein descent','x1','x2');
drawnow;

NbIter = length(x_history);
printf('Number of call to the objective function = %d\n',NbIter);

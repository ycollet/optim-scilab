//functoplot='rosenbrock';
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
functoplot='abspow';
//functoplot='michalewicz';
//functoplot='levy5';
//functoplot='levy13';

deff('y=fobj(x)','y='+functoplot+'(x)');

// example of use of the simulated tempering method

Proba_start = [0.1 0.3 0.6 0.9];
ItMax    = 1000;
It_Pre   = 100;
st_proba = 0.05;
NbTemp   = 20;

//////////////////////////////////////////

x_min = eval('get_min_bound_'+functoplot+'()');
x_max = eval('get_max_bound_'+functoplot+'()');

x0 = (x_max - x_min) .* rand(size(x_min,1),size(x_min,2)) + x_min;

x = x_min(1):(x_max(1) - x_min(1))/ 10:x_max(1);
y = x_min(2):(x_max(2) - x_min(2))/ 10:x_max(2);

[X,Y]=meshgrid(x,y);
for i=1:size(X,1)
  for j=1:size(X,2)
    Z(i,j) = eval(functoplot+'([X(i,j) Y(i,j)])');
  end
end

/////////////////////////
// Simulated Tempering //
/////////////////////////

printf('ST: geometrical decrease temperature law\n');

delta = [-0.1*(x_max - x_min)' 0.1*(x_max - x_min)']';

T0 = compute_initial_temp(x0, fobj, Proba_start, It_Pre, neigh_func_default, neigh_func_param = delta);
printf('Initial temperatore T0 = '); disp(T0);
T = T0;

//// Generation of the temperature list
//T(NbTemp) = T0;
//for i=NbTemp-1:-1:1
//  T(i) = T(i+1) / 2.0;
//end

[x_opt, temp_list, x_history] = optim_st(x0, fobj, ItMax, T, st_proba, Log = %T, para_neigh_func = delta);

printf('Solution found:'); disp(x_opt');
printf('value of the objective function = %f\n', fobj(x_opt));
printf('Optimal solution:'); disp(eval('get_opti_'+functoplot+'()')');

drawlater;

xset('fpf',' ');
contour(x,y,Z,10);
_axes = get("current_axes");
_axes.data_bounds = [x_min(1) x_max(1) x_min(2) x_max(2)];

for i=1:length(x_history)-1
  plot([x_history(i)(1) x_history(i+1)(1)],[x_history(i)(2) x_history(i+1)(2)],'k-');
  plot(x_history(i)(1),x_history(i)(2),'ro');
end
plot(x_history($)(1),x_history($)(2),'ro');
drawnow;

scf();
xtitle('Temperature evolution','Iteration','Temperature');
t = 1:length(temp_list);
plot(t, temp_list,'k-');


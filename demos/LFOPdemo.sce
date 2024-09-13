FuncName = 'rosenbrock';

function y = df_rosenbrock(x)
y(1) = - 100*4*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
y(2) = 100*2*(x(2)-x(1)^2);
endfunction

deff('y=df_rosenbrock_estim(x)','y=derivative(rosenbrock,x)''');

x0 = [-1.2; 1];
// Parameters for try 1
delta_t_1    = 0.5;  // the temporal step size for the resolution of the dynamical system
delta_step_1 = 0.1;  // maximal allowable step size
GradTOL_1    = 1e-5; // we stop the algorithm if the gradient is below the GradTOL level.
delta_inc_1  = 0.01; // increase factor of the dilatation coeff of delta_t
m_1          = 3;    // id ak+1'.ak is non positive during m iterations then delta_t = delta_t / 2
                     // and we restart from (xk+xk+1) / 2
p_start_1    = 1.1;  // starting value of the dilatation coefficient: 1.01 by default

// Parameters for try 2
delta_t_2    = 0.5;  // the temporal step size for the resolution of the dynamical system
delta_step_2 = 0.1;  // maximal allowable step size
GradTOL_2    = 1e-5; // we stop the algorithm if the gradient is below the GradTOL level.
delta_inc_2  = 0.01; // increase factor of the dilatation coeff of delta_t
m_2          = 3;    // id ak+1'.ak is non positive during m iterations then delta_t = delta_t / 2
                     // and we restart from (xk+xk+1) / 2
p_start_2    = 2.1;  // starting value of the dilatation coefficient: 1.01 by default


// delta_t    : a une influence au début de la recherche. Il vaut mieux le choisir petit
// delta_step : à choisir en fonction de la taille de l'espace de recherche : ~(max(x(i)) - min(x(i)))/100 semble être un bon départ. Si le pas max
//              est trop grand,on remplace le pas courant par un pas beaucoup  trop grand en fin de recherche.
// GradTOL    : si trop grand alors arrêt prématuré (1e-5 / 1e-1)
// delta_inc  : si trop grand, la taille de delta_t remonte trop vite vers la fin de la recherche.
// p_start    : si trop grand, la taille de delta_t est réinitialisée périodiquement à une valeur trop grand au cours de la recherche.
// m          : si trop grand, pas assez d'adaptation du pas.
MaxEvalFunc = 400;
step_x_history = 1;
Plot = %T;
offset = 0.25;
Min = get_min_bound_rosenbrock();
Max = get_max_bound_rosenbrock();

if (Plot) then
  clf;
  x = Min(1):(Max(1) - Min(1))/20:Max(1);
  y = Min(2):(Max(2)-Min(2))/20:Max(2);
  for i=1:size(x,2)
    for j=1:size(y,2)
      Z(i,j) = rosenbrock([x(i) y(j)]);
    end
  end
  xset('fpf',' ');
  contour(x,y,Z, 10);
end

[x_opt, x_history] = optim_LFOP(df_rosenbrock, x0, delta_t_1, delta_step_1, m_1, delta_inc_1, GradTOL_1, MaxEvalFunc, %T, p_start_1);

printf('x_opt = '); disp(x_opt);
printf('f(x_opt) = %f\n', evstr(FuncName+'(x_opt)'));

if (Plot) then
  wId = waitbar(0,'Drawing results');
  drawlater;
  for i=1:step_x_history:length(x_history)-1
    if (modulo(i, ceil((length(x_history)-1) / (100/step_x_history)))==0) then
      waitbar(floor(1000*i/length(x_history))/1000,wId);
    end
    plot([x_history(i)(1) x_history(i+1)(1)], [x_history(i)(2) x_history(i+1)(2)], 'k-');
    plot(x_history(i+1)(1), x_history(i+1)(2), 'ro');
  end
  xtitle('LFOP - comparison of parameter influence','x1','x2');
  drawnow;
  winclose(wId);
end

[x_opt, x_history] = optim_LFOP(df_rosenbrock_estim, x0, delta_t_2, delta_step_2, m_2, delta_inc_2, GradTOL_2, MaxEvalFunc, %T, p_start_2);

printf('x_opt = '); disp(x_opt);
printf('f(x_opt) = %f\n', evstr(FuncName+'(x_opt)'));

if (Plot) then
  wId = waitbar(0,'Drawing results');
  drawlater;
  for i=1:step_x_history:length(x_history)-1
    if (modulo(i, ceil((length(x_history)-1) / (100/step_x_history)))==0) then
      waitbar(floor(1000*i/length(x_history))/1000,wId);
    end
    plot([x_history(i)(1) x_history(i+1)(1)], [x_history(i)(2)-offset x_history(i+1)(2)-offset], 'k-.');
    plot(x_history(i+1)(1), x_history(i+1)(2)-offset, 'go');
  end
  drawnow;
  winclose(wId);
  legends(['Try 1','Try 2'],[3,5],1);
end


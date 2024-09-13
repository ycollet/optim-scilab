//FuncName = 'rosenbrock';
//FuncName = 'mccormic';
FuncName = 'sixhumpcamelback';

deff('y=f(x)','y='+FuncName+'(x)');
deff('y=df(x)','y=derivative('+FuncName+',x)''');
deff('y=get_min_bound()','y=get_min_bound_'+FuncName+'()');
deff('y=get_max_bound()','y=get_max_bound_'+FuncName+'()');

Min = get_min_bound();
Max = get_max_bound();

x0 = [-1.2; 1];

ItMX = 400;
XTol = 1e-6;
GradTol = 1e-6;
rho  = min(Max-Min)/100;
Plot = %T;
Log  = %T;



if (Plot) then
  clf;
  x = Min(1):(Max(1) - Min(1))/20:Max(1);
  y = Min(2):(Max(2) - Min(2))/20:Max(2);
  for i=1:size(x,2)
    for j=1:size(y,2)
      Z(i,j) = rosenbrock([x(i) y(j)]);
    end
  end
  xset('fpf',' ');
  contour(x,y,Z, 10);
end

[x_opt, x_history] = optim_sqsd(f, df, x0, ItMX, XTol, GradTol, rho, Log);

printf('x_opt = '); disp(x_opt);
printf('f(x_opt) = %f\n', evstr(FuncName+'(x_opt)'));

if (Plot) then
  wId = waitbar(0,'Drawing results');
  drawlater;
  for i=1:length(x_history)-1
    if (modulo(i, ceil((length(x_history)-1) / 100))==0) then
      waitbar(floor(1000*i/length(x_history))/1000,wId);
    end
    plot([x_history(i)(1) x_history(i+1)(1)], [x_history(i)(2) x_history(i+1)(2)], 'k-');
    plot(x_history(i+1)(1), x_history(i+1)(2), 'ro');
  end
  xtitle('SQSD - Function '+FuncName,'x1','x2');
  drawnow;
  winclose(wId);
end


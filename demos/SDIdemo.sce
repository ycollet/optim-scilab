FuncName1 = 'treccani';
FuncName2 = 'matyas';

deff('[f1, f2] = func(x)','f1 = '+FuncName1+'(x); ...
                           f2 = '+FuncName2+'(x);');
                           
x0 = [-1.2; 1];
noise_type  = ['uniform';'uniform'];
noise_param = [0.4, 0.4; -0.2, -0.2];
target   = [0, 0];
weight   = [1, 1];
NbLoops  = 10;
NbPoints = 10;
Plot     = %T;
Min = [-1.5, -1.5];
Max = [1.5, 1.5];

[InitPoint1, InitPoint2] = func(x0); disp(InitPoint1);
printf('Initial point: %f %f\n', InitPoint1, InitPoint2);
printf('x0 ='); disp(x0);

[x_opt, x_history, x_pool_history] = optim_sdi(func, x0, noise_type, noise_param, target, weight, NbLoops, NbPoints);

[FinalPoint1, FinalPoint2] = func(x_opt); disp(FinalPoint1);
printf('Final point: %f %f\n', FinalPoint1, FinalPoint2);
printf('x_opt = '); disp(x_opt);
printf('Improvement: %f \% %f \%\n', (InitPoint1 - FinalPoint1)./InitPoint1 * 100, (InitPoint2 - FinalPoint2)./InitPoint2 * 100);

Min(1) = min([Min(1) x_opt(1)]);
Min(2) = min([Min(2) x_opt(2)]);
Max(1) = max([Max(1) x_opt(1)]);
Max(2) = max([Max(2) x_opt(2)]);

if (Plot) then
  clf();
  drawlater;
  x = Min(1):(Max(1) - Min(1))/20:Max(1);
  y = Min(2):(Max(2) - Min(2))/20:Max(2);
  for i=1:size(x,2)
    for j=1:size(y,2)
      Z1(i,j) = evstr(FuncName1+'([x(i) y(j)])');
      Z2(i,j) = evstr(FuncName2+'([x(i) y(j)])');
    end
  end

  subplot(2,1,1);
  contour(x, y, Z1, 10);
  xset('fpf',' ');
  xtitle('SDI: f1 = '+FuncName1+' - target = 0 - weight = 1','x1','x2');
  
  subplot(2,1,2);
  contour(x, y, Z2, 10);
  xset('fpf',' ');
  xtitle('SDI: f2 = '+FuncName2+' - target = 0 - weight = 1','x1','x2');
  
  wId = waitbar(0,'Drawing results');
  
  subplot(2,1,1);
  // Plot the pools
  for i=1:length(x_pool_history)
    waitbar(i/length(x_pool_history), wId);
    for j=1:length(x_pool_history(i))
      plot(x_pool_history(i)(j)(1), x_pool_history(i)(j)(2), 'go');
      handle_pt = gce();
      handle_pt.children.data(3) = eval(FuncName1+'([x_pool_history(i)(j)(1), x_pool_history(i)(j)(2)])');
    end
  end

  subplot(2,1,2);
  // Plot the pools
  for i=1:length(x_pool_history)
    waitbar(i/length(x_pool_history), wId);
    for j=1:length(x_pool_history(i))
      plot(x_pool_history(i)(j)(1), x_pool_history(i)(j)(2), 'go');
      handle_pt = gce();
      handle_pt.children.data(3) = eval(FuncName2+'([x_pool_history(i)(j)(1), x_pool_history(i)(j)(2)])');
    end
  end

  drawnow;
  winclose(wId);
end


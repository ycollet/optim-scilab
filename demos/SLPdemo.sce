//deff('y=f(x)','y=sum(x.^2)');
//deff('y=df(x)','y=2*x');
//deff('y=ineqconstraint(x)','y(1)=-[1 2]*x-[3 4]*x.^2+1; y(2) = -[1 4]*x-[4 5]*x.^2+2');
//deff('y=df_ineqconstraint(x)','y=-[1 1;2 4]');
//eqconstraint    = [];
//df_eqconstraint = [];
//upper       = [2;2];
//lower       = [-2; -2];
//x0          = [2; 2];

//deff('y=f(x)','y=x(1)^2 + x(2)^2 - 16*x(1) - 10*x(2)');
//deff('y=df(x)','y=[2*x(1) - 16 ; 2*x(2) - 10]');
//deff('y=ineqconstraint(x)','y(1) = - 11 + x(1)^2 - 6*x(1) + 4*x(2); ...
//                        y(2) = - x(1)*x(2) + 3*x(2) + exp(x(1) - 3) - 1; ...
//                        y(3) = - x(1); ...
//                        y(4) = - x(2);');
//deff('y=df_ineqconstraint(x)','y(1,1) = 2*x(1) - 6; ...
//                           y(2,1) = 4; ...
//                           y(1,2) = -x(2) + exp(x(1) - 3); ...
//                           y(2,2) = -x(1) + 3; ...
//                           y(1,3) = -1; ...
//                           y(2,3) = 0; ...
//                           y(1,4) = 0; ...
//                           y(2,4) = -1;');
//eqconstraint    = [];
//df_eqconstraint = [];
//upper       = [15;9];
//lower       = [0; 0];
//x0          = [4; 3];

//deff('y=f(x)','y=5*x(1)^2-3*x(2)^2');
//deff('y=df(x)','y(1) = 10*x(1); ...
//                y(2) = -6*x(2);');
//deff('y=ineqconstraint(x)','y(1) = -x(1); ...
//                        y(2) = -x(2);');
//deff('y=df_ineqconstraint(x)','y(1,1)=-1; ...
//                           y(2,1)=0; ...
//                           y(1,2)=0; ...
//                           y(2,2)=-1;');
//eqconstraint    = [];
//df_eqconstraint = [];
//upper = [4;4];
//lower = [-4;-4];
//x0 = [1;1];

//deff('y=f(x)','y=x(1)*x(2)');
//deff('y=df(x)','y(1) = x(2); ...
//                y(2) = x(1);');
//deff('y=ineqconstraint(x)','y(1) = -25 + x(1)^2 + x(2)^2;');
//deff('y=df_ineqconstraint(x)','y(1,1) = 2*x(1); ...
//                           y(2,1) = 2*x(2);');
//eqconstraint    = [];
//df_eqconstraint = [];
//upper = [4;4];
//lower = [-4;-4];
//x0    = [0.1;0.1]; // Feasible starting point
////x0    = [4;4]; // Non feasible starting point

deff('y=f(x)','y=4*x(1)- x(2)^2 - 12');
deff('y=df(x)','y(1) = 4; ...
                y(2) = -2*x(2);');
deff('y=eqconstraint(x)','y(1) = 20 - x(1)^2 - x(2)^2');
deff('y=df_eqconstraint(x)','y(1,1) = -2*x(1); ...
                             y(2,1) = -2*x(2);');
deff('y=ineqconstraint(x)','y(1) = - 10*x(1) + x(1)^2 - 10*x(2) + x(2)^2 + 34; ...
                            y(2) = - x(1); ...
                            y(3) = - x(2);');
deff('y=df_ineqconstraint(x)','y(1,1) = -10 + 2*x(1); ...
                               y(2,1) = -10 + 2*x(2); ...
                               y(1,2) = -1; ...
                               y(2,2) = 0; ...
                               y(1,3) = 0; ...
                               y(2,3) = -1;');
upper = [14;14];
lower = [-14;-14];
//x0    = [2;4]; // Feasible starting point
x0    = [12;4]; // non-feasible starting point

//deff('y=f(x)','y=x(1)^2+x(2)^2');
//deff('y=df(x)','y(1) = 2*x(1); ...
//                y(2) = 2*x(2);');
//deff('y=eqconstraint(x)','y(1) = - x(1)^2 - x(2)^2 + 9*x(2) - 4.25');
//deff('y=df_eqconstraint(x)','y(1,1) = -2*x(1); ...
//                             y(2,1) = -2*x(2) + 9;');
//ineqconstraint = [];
//df_ineqconstraint = [];
//upper = [4;4];
//lower = [-4;-4];
//x0    = [2;3.9]; // Non feasible starting point

//// Test with some parameters for function f
//deff('y=f(x,param1,param2)','y=sum(x.^2); ...
//                             if ~isdef(''param1'',''local'') then param1 = -1; end; ...
//                             if ~isdef(''param2'',''local'') then param2 = -2; end; ...
//                             printf(''param1 = %d\n'',param1); ...
//                             printf(''param2 = %d\n'',param2);');

//deff('y=f(x)','y=sum(x.^2)');
//deff('y=df(x)','y=2*x');
//deff('y=ineqconstraint(x)','y(1)=-[1 2]*x-[3 4]*x.^2+1; y(2) = -[1 4]*x-[4 5]*x.^2+2');
//deff('y=df_ineqconstraint(x)','y(1,1) = -1 - 6*x(1); ...
//                               y(1,2) = -2 - 8*x(2); ...
//                               y(2,1) = -1 - 8*x(1); ...
//                               y(2,2) = -4 - 10*x(2);');
//eqconstraint    = [];
//df_eqconstraint = [];
//upper       = [2;2];
//lower       = [-2; -2];
//x0          = [2; 2];

ItMX        = 1000;
MaxEvalFunc = 10000;
StepPlot    = 1;
//delta_ml    = 0.005;
delta_ml    = []; // To test the automatic sizing of the move limits
Plot        = %T;
DispNum     = %T;
Log         = %T;
MaxMinStepCount = 5;
//Weight = [1 1 1 1];
Weight = [];
XTOL = 1e-4;
ITOL = 1e-4;
ETOL = 1e-4;

//if (Plot & max(size(x0))==2)     then Plot2 = %T;
//elseif (Plot & max(size(x0))==1) then Plot1 = %T;
//else
//  Plot1 = %F;
//  Plot2 = %F;
//end
Plot2 = %T; Plot1 = %F;

[x_opt, x_history, ml_history] = optim_slp(f, df, ...
                                           ineqconstraint, df_ineqconstraint, ...
                                           eqconstraint, df_eqconstraint, ...
                                           x0, ItMX, MaxEvalFunc, delta_ml, upper, lower, Weight, MaxMinStepCount, Log, XTOL, ITOL, ETOL);

printf('Initial point\n');
if (ineqconstraint~=[]) then
  printf('Value of the inequality constraints (must be all negatives or null) :'); disp(ineqconstraint(x0)');
end
if (eqconstraint~=[]) then
  printf('Value of the equality constraints (must be all null) :'); disp(eqconstraint(x0)');
end
printf('Value of the objective function : %f\n', f(x0));

printf('Solution:\n');
if (ineqconstraint~=[]) then
  printf('Value of the inequality constraints (must be all negatives or null) :'); disp(ineqconstraint(x_opt)');
end
if (eqconstraint~=[]) then
  printf('Value of the equality constraints (must be all null) :'); disp(eqconstraint(x_opt)');
end
printf('Value of objective function = %f\n', f(x_opt));
printf('Solution found:'); disp(x_opt');

if (Plot) then
  scf();
  T = 1:length(x_history);
  ListF = [];
  for i=1:length(x_history);
    ListF(i) = f(x_history(i));
  end
  // We plot the objective function
  plot(T,ListF,'k');
  xtitle('Evolution of the objective function','Iteration','F');
  // We plot the equality and inequality constraints function
  if (eqconstraint~=[]) then
    scf();
    ListEC = [];
    for i=1:length(x_history)
      ListEC(:,i) = eqconstraint(x_history(i));
    end
    for i=1:size(ListEC,1)
      plot(1:length(x_history),ListEC(i,:),'k');
    end
    xtitle('Evolution of the equality constraints functions','Iteration','H');
  end
  if (ineqconstraint~=[]) then
    scf();
    ListIC = [];
    for i=1:length(x_history)
      ListIC(:,i) = ineqconstraint(x_history(i));
    end
    for i=1:size(ListIC,1)
      plot(1:length(x_history),ListIC(i,:),'k');
    end
    xtitle('Evolution of the inequality constraints functions','Iteration','H');
  end
end

if (Plot1) then
  scf();
  drawlater;
  x = lower(1):(upper(1) - lower(1))/20:upper(1);
  for i=1:length(x)
    Z(i) = f(x(i));
  end
  plot(x,Z);

  wId = waitbar(0,'Drawing results');
  for i=1:StepPlot:length(x_history)
    if (modulo(i/StepPlot, ceil((length(x_history)/StepPlot) / 10))==0) then
      waitbar(floor(1000*i/length(x_history))/1000,wId);
    end
    plot(x_history(i)(1), f(x_history(i)(1)), 'ro');
    if (i~=length(x_history)) then
      FrameColor = 'g-';
    else
      FrameColor = 'b-';
    end
    plot([x_history(i)(1) + ml_history(i)(1) x_history(i)(1) - ml_history(i)(1)],[f(x_history(i)(1)) f(x_history(i)(1))], FrameColor);
    if (DispNum) then
      xstring(x_history(i)(1) - ml_history(i)(1),f(x_history(i)(1)), string(i));
    end
  end
  xtitle('SLP','x');
  legends(['Move limits','Solution found'],[3,5],1);
  drawnow;
  winclose(wId);
end

if (Plot2) then
  scf();
  drawlater;
  x = lower(1):(upper(1) - lower(1))/20:upper(1);
  y = lower(2):(upper(2) - lower(2))/20:upper(2);
  for i=1:size(x,2)
    for j=1:size(y,2)
      Z(i,j) = f([x(i) y(j)]);
    end
  end
  xset('fpf',' ');
  contour(x,y,Z, 10);
  for i=1:StepPlot:length(x_history)
    Rect(1,i) = x_history(i)(1) - ml_history(i)(1);
    Rect(2,i) = x_history(i)(2) + ml_history(i)(2);
    Rect(3,i) = 2*ml_history(i)(1);
    Rect(4,i) = 2*ml_history(i)(2);
    X(i) = x_history(i)(1);
    Y(i) = x_history(i)(2);
    Color(i) = -3;
    if i==length(x_history) then Color(i) = -11; end
    if (DispNum) then
      xstring(x_history(i)(1) - ml_history(i)(1), x_history(i)(2) - ml_history(i)(2), string(i));
    end
  end
  xrects(Rect,Color);
  plot(X, Y, 'ro');
  xtitle('SLP','x1','x2');
  legends(['Move limits','Solution found'],[3,5],1);
  drawnow;
end


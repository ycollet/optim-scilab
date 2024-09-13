lines(0);

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
//upper = [4;4];
//lower = [-4;-4];
//x0 = [1;1];

//deff('y=f(x)','y=x(1)*x(2)');
//deff('y=df(x)','y(1) = x(2); ...
//                y(2) = x(1);');
//deff('y=ineqconstraint(x)','y(1) = -25 + x(1)^2 + x(2)^2;');
//deff('y=df_ineqconstraint(x)','y(1,1) = 2*x(1); ...
//                           y(2,1) = 2*x(2);');
//upper = [4;4];
//lower = [-4;-4];
////x0    = [0.1;0.1]; // Feasible starting point
//x0    = [4;4]; // Non feasible starting point

deff('y=f(x)','y=sum(x.^2)');
deff('y=df(x)','y=2*x');
deff('y=ineqconstraint(x)','y(1)=-[1 2]*x-[3 4]*x.^2+1; y(2) = -[1 4]*x-[4 5]*x.^2+2');
deff('y=df_ineqconstraint(x)','y(1,1) = -1 - 6*x(1); ...
                               y(1,2) = -2 - 8*x(2); ...
                               y(2,1) = -1 - 8*x(1); ...
                               y(2,2) = -4 - 10*x(2);');
eqconstraint    = [];
df_eqconstraint = [];
upper       = [2;2];
lower       = [-2; -2];
x0          = [2; 2];

ItMX            = 100;
ls_ItMX         = 100;
MaxEvalFunc     = 100;
MaxMinStepCount = 10;
h        = 0.00125;
SubTOL   = 1.0e-2;
Log      = %F;
Plot     = %T;
DispNum  = %F;
StepPlot = 1;
epsilon  = -0.01;
Theta_0  = 1.0;


[x_opt, x_history] = optim_feasdir(f, df, ineqconstraint, df_ineqconstraint, ...
                                   x0, ItMX, MaxEvalFunc, upper, lower, Theta_0, epsilon, h, MaxMinStepCount, ...
                                   SubTOL, ls_ItMX, Log);

printf('Initial point\n');
printf('Value of the inequality constraints (must be all negatives or null) :'); disp(ineqconstraint(x0)');
printf('Value of the objective function : %f\n', f(x0));

printf('Solution:\n');
printf('Value of the inequality constraints (must be all negatives or null) :'); disp(ineqconstraint(x_opt)');
printf('Value of objective function = %f\n', f(x_opt));
printf('Solution found:'); disp(x_opt');

if (Plot) then
  scf();
  EvalFunc  = 0;
  IndexList = [];
  ListF     = [];
  Index     = 1;
  for i=1:length(x_history);
    IndexList($+1) = Index;
    disp(x_history(i)($)');
    printf('f(x_history(i)($)) = %f\n', f(x_history(i)($)));
    printf('g(x_history(i)($)) = '); disp(ineqconstraint(x_history(i)($))');
    EvalFunc = EvalFunc + length(x_history(i));
    for j=1:length(x_history(i))
      ListF(Index) = f(x_history(i)(j));
      Index = Index + 1;
    end
  end
  T = 1:length(ListF);
  // We plot the objective function
  plot(T,ListF,'k');
  plot(IndexList, ListF(IndexList),'ro');
  xtitle('Evolution of the objective function','Iteration','F');
  // We plot the equality and inequality constraints function
  if (ineqconstraint~=[]) then
    scf();
    IndexList = [];
    ListIC    = [];
    Index     = 1;
    for i=1:length(x_history);
      IndexList($+1) = Index;
      EvalFunc = EvalFunc + length(x_history(i));
      for j=1:length(x_history(i))
        ListIC(:,Index) = ineqconstraint(x_history(i)(j));
        Index = Index + 1;
      end
    end
    T = 1:size(ListIC,2);
    for i=1:size(ListIC,1)
      plot(T,ListIC(i,:),'k');
      plot(IndexList,ListIC(i,IndexList),'ro');
    end
    xtitle('Evolution of the inequality constraints functions','Iteration','H');
  end
  printf('Number of evaluation of the objective function: %d\n', EvalFunc);
end


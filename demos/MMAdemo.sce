lines(0);

deff('y=f(x)','y=x(1)^2 + x(2)^2 - 16*x(1) - 10*x(2)');
deff('y=df(x)','y=[2*x(1) - 16 ; 2*x(2) - 10]');
deff('y=ineqconstraint(x)','y(1) = - 11 + x(1)^2 - 6*x(1) + 4*x(2); ...
                        y(2) = - x(1)*x(2) + 3*x(2) + exp(x(1) - 3) - 1; ...
                        y(3) = - x(1); ...
                        y(4) = - x(2);');
deff('y=df_ineqconstraint(x)','y(1,1) = 2*x(1) - 6; ...
                           y(2,1) = 4; ...
                           y(1,2) = -x(2) + exp(x(1) - 3); ...
                           y(2,2) = -x(1) + 3; ...
                           y(1,3) = -1; ...
                           y(2,3) = 0; ...
                           y(1,4) = 0; ...
                           y(2,4) = -1;');
eqconstraint    = [];
df_eqconstraint = [];
upper       = [15;9];
lower       = [0; 0];
x0          = [4; 3];

//deff('y=f(x)','y=4*x(1)- x(2)^2 - 12');
//deff('y=df(x)','y(1) = 4; ...
//                y(2) = -2*x(2);');
//deff('y=eqconstraint(x)','y(1) = 25 - x(1)^2 - x(2)^2');
//deff('y=df_eqconstraint(x)','y(1,1) = -2*x(1); ...
//                             y(2,1) = -2*x(2);');
//deff('y=ineqconstraint(x)','y(1) = - 10*x(1) + x(1)^2 - 10*x(2) + x(2)^2 + 34; ...
//                            y(2) = - x(1); ...
//                            y(3) = - x(2);');
//deff('y=df_ineqconstraint(x)','y(1,1) = -10 + 2*x(1); ...
//                               y(2,1) = -10 + 2*x(2); ...
//                               y(1,2) = -1; ...
//                               y(2,2) = 0; ...
//                               y(1,3) = 0; ...
//                               y(2,3) = -1;');
//
//upper = [14;14];
//lower = [-14;-14];
//x0    = [12;12]; // Feasible starting point
////x0    = [-12;12]; // Feasible starting point

ItMX        = 25;
MaxEvalFunc = 1000;
Log         = %T;
XTol        = 1e-6;
DispNum     = %F;

[x_opt, x_history] = optim_mma(f, df, ineqconstraint, df_ineqconstraint, eqconstraint, df_eqconstraint, x0, upper, lower, ItMX, MaxEvalFunc, delta_step, Log);

printf('\n----------\n\n');

printf('Solution for each penalty coefficient:\n');
for i=1:length(x_history)
  printf('Step %d: optimal parameters:', i); disp(x_history(i)')
  printf('Step %d: objective function value = %f\n', i, f(x_history(i))); 
  printf('Step %d: constraint values :', i); disp(ineqconstraint(x_history(i))');  
end
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
for i=1:length(x_history)
  plot(x_history(i)(1), x_history(i)(2), 'ro');
  if (DispNum) then
    xstring(x_history(i)(1), x_history(i)(2), string(i));
  end
end
plot(x_history($)(1), x_history($)(2), 'gx');
xtitle('MMA','x1','x2');
legends(['Solution found','Intermediate solutions'],[3,5],1);
drawnow;

printf('Solution:\n');
printf('optimal parameters:'); disp(x_opt')
printf('objective function value = %f\n', f(x_opt)); 
printf('constraint values :'); disp(ineqconstraint(x_opt)');

printf('\n----------\n\n');

printf('Solution for the initial solution:\n');
printf('optimal parameters:'); disp(x0')
printf('objective function value = %f\n', f(x0)); 
printf('constraint values :'); disp(ineqconstraint(x0)');




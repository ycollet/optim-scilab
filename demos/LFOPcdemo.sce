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
//upper = [14;14];
//lower = [-14;-14];
//x0    = [12;12]; // Feasible starting point

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

MaxEvalFunc = 1000;
delta_t     = 0.5;
delta_step  = 0.1;
m           = 3;
delta_inc   = 0.01;
GradTOL     = 1e-4;
p_start     = 1.1;
Log         = %T;


[x_opt, x_history] = optim_LFOPc(df, ineqconstraint, df_ineqconstraint, eqconstraint, df_eqconstraint, ...
                                 x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);

printf('Initial point:'); disp(x0');
if (ineqconstraint~=[]) then
  printf('inequality constraints:'); disp(ineqconstraint(x0)');
end
if (eqconstraint~=[]) then
  printf('equality constraints:'); disp(eqconstraint(x0)');
end
printf('Objective function value = %f\n', f(x0));

printf('Solution:'); disp(x_opt');
if (ineqconstraint~=[]) then
  printf('inequality constraints:'); disp(ineqconstraint(x_history($-2))');
end
if (eqconstraint~=[]) then
  printf('equality constraints:'); disp(eqconstraint(x_history($-2))');
end
printf('Objective function value = %f\n', f(x_history($-2)));

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

plot(x_history(1)(1), x_history(1)(2), 'ro');
for i=1:length(x_history)
  plot2d(x_history(i)(1), x_history(i)(2), style = 0);
end
plot(x_history($-2)(1), x_history($-2)(2), 'go');
xtitle('LFOPc demo','x1','x2');
drawnow;

//deff('y=f(x)','y=sum(x.^2)');
//deff('y=df(x)','y=2*x');
//deff('y=ineqconstraint(x)','y(1)=-[1 2]*x-[3 4]*x.^2+1; y(2) = -[1 4]*x-[4 5]*x.^2+2');
//deff('y=df_ineqconstraint(x)','y(1,1) = -1 - 6*x(1); ...
//                               y(1,2) = -2 - 8*x(2); ...
//                               y(2,1) = -1 - 8*x(1); ...
//                               y(2,2) = -4 - 10*x(2);y=y'';');
//eqconstraint    = [];
//df_eqconstraint = [];
//upper       = [2;2];
//lower       = [-2; -2];
//x0          = [2; 2];

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

ItMX    = 10000;
NbLoop  = 3;
Log     = %T;
DispNum = %T;
r_p     = 5;
eta     = 5;
r_max   = 1000;

[x_opt, x_history] = optim_laau(f, df, ...
                                ineqconstraint, df_ineqconstraint, ...
                                eqconstraint, df_eqconstraint, x0, ItMX, NbLoop, Log, r_p, eta, r_max);

printf('Solution:\n');
printf('optimal parameters:'); disp(x_opt')
printf('objective function value = %f\n', f(x_opt)); 
printf('constraint values :'); disp(ineqconstraint(x_opt)');

printf('\n----------\n\n');

printf('Solution for the initial solution:\n');
printf('optimal parameters:'); disp(x0')
printf('objective function value = %f\n', f(x0)); 
printf('constraint values :'); disp(ineqconstraint(x0)');

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
  X(i) = x_history(i)(1);
  Y(i) = x_history(i)(2);
  Color(i) = -3;
  if i==length(x_history) then Color(i) = -11; end
  if (DispNum) then
    xstring(x_history(i)(1), x_history(i)(2), string(i));
  end
end
plot(X, Y, 'ro');
xtitle('SUMT','x1','x2');
legends(['Solution found'],[3],1);
drawnow;


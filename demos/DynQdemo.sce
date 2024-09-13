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

////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////

//deff('y=f(x)','y=x(1)^2+x(2)^2');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=eqconstraint(x)','y(1) = - x(1)^2 - x(2)^2 + 9*x(2) - 4.25');
//deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');
//ineqconstraint = [];
//df_ineqconstraint = [];
//upper = [4;4];
//lower = [-4;-4];
//x0    = [2;3.9]; // Non feasible starting point

////////////////////////////////////////////////////////////////////////

//// constr_pb_2
//deff('y=f(x)','y=x(1)^2 + x(2)^2 - 16*x(1) - 10*x(2)');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=ineqconstraint(x)','y(1) = - 11 + x(1)^2 - 6*x(1) + 4*x(2); ...
//                        y(2) = - x(1)*x(2) + 3*x(2) + exp(x(1) - 3) - 1; ...
//                        y(3) = - x(1); ...
//                        y(4) = - x(2);');
//deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');
//eqconstraint = [];
//df_eqconstraint = [];
//
//upper       = [15;9];
//lower       = [0; 0];
//x0          = [4; 3];

////////////////////////////////////////////////////////////////////////
//constr_pb_1
deff('y=f(x)','y=4*x(1)- x(2)^2 - 12');
deff('y=df(x)','y = derivative(f,x)''');
deff('y=ineqconstraint(x)','y(1,1) = - 10*x(1) + x(1)^2 - 10*x(2) + x(2)^2 + 34; ...
                            y(2,1) = - x(1); ...
                            y(3,1) = - x(2);');
deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');
deff('y=eqconstraint(x)','y(1) = 20 - x(1)^2 - x(2)^2');
deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');

upper = [14;14];
lower = [-14;-14];
x0    = [-12;-12]; // Feasible starting point
//x0    = [1.5;1.5]; // Feasible starting point

////////////////////////////////////////////////////////////////////////

//deff('y=f(x)','y=sum(x.^2)');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=ineqconstraint(x)','y(1)=-[1 2]*x-[3 4]*x.^2+1; y(2) = -[1 4]*x-[4 5]*x.^2+2');
//deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');
//eqconstraint    = [];
//df_eqconstraint = [];
//upper       = [2;2];
//lower       = [-2; -2];
//x0          = [1.99; 1.99];
////x0          = [0; 0];

////////////////////////////////////////////////////////////////////////

////funcname = 'G1';
////funcname = 'G2';
////funcname = 'G3';
////funcname = 'G4';
////funcname = 'G5';
////funcname = 'G6';
////funcname = 'G7';
////funcname = 'G8';
////funcname = 'G9';
////funcname = 'G10';
////funcname = 'G11';
////funcname = 'G12';
////funcname = 'G13';
////funcname = 'weldedbeam';
//funcname = 'pressurevessel';
////funcname = 'tensioncompr'; // ??
////funcname = 'elec_on_sph'; // ??
////funcname = 'larg_sm_poly'; // ??
////funcname = 'vv_open'; // ??
////funcname = 'himmelblau_1';
////funcname = 'himmelblau_3';
////funcname = 'himmelblau_4';
////funcname = 'himmelblau_4a';
////funcname = 'himmelblau_5';
////funcname = 'himmelblau_8';
////funcname = 'himmelblau_9';
////funcname = 'himmelblau_10';
////funcname = 'himmelblau_11';
////funcname = 'himmelblau_12';
////funcname = 'himmelblau_13';
////funcname = 'himmelblau_14';
////funcname = 'himmelblau_15';
////funcname = 'himmelblau_16';
////funcname = 'himmelblau_24';
//
//x0 = eval(funcname+'_x_init()');
//upper = eval('get_min_'+funcname+'()');
//lower = eval('get_max_'+funcname+'()');
//deff('y=f(x)','y = '+funcname+'_obj(x)');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=ineqconstraint(x)','y = '+funcname+'_ineq_constr(x)');
//deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');
//deff('y=eqconstraint(x)','y = '+funcname+'_eq_constr(x)');
//deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');

////////////////////////////////////////////////////////////////////////

////function [value] = elec_on_sph_obj(x)
////function [diff_bound] = elec_on_sph_constr(x)
//N = 10;
//upper = ones(3*N,1);
//lower = -ones(3*N,1);
//x0 = (upper - lower) .* rand(size(upper,1),size(upper,2)) + lower;
//
//deff('y=f(x)','y = elec_on_sph_obj(matrix(x,N,3))');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=ineqconstraint(x)','y = elec_on_sph_ineq_constr(matrix(x,N,3))');
//deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');

////////////////////////////////////////////////////////////////////////

////function [value] = larg_sm_poly_obj(r, theta)
////function [diff_bound] = larg_sm_poly_constr(r, theta)
//// 20 variables r, 20 variables theta
//N = 6;
//upper_r = ones(N,1);
//lower_r = zeros(N,1);
//upper_theta = %pi*ones(N,1);
//lower_theta = zeros(N,1);
//upper = [upper_r' upper_theta']';
//lower = [lower_r' lower_theta']';
//r0 = 0.9*ones(N,1);
//theta0 = (0:%pi/(N-1):%pi)';
//x0 = [r0' theta0']';
//
//deff('y=f(x)','y = larg_sm_poly_obj(x(1:N),x(N+1:$))');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=ineqconstraint(x)','y = larg_sm_poly_ineq_constr(x(1:N),x(N+1:$))');
//deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');

/////////////////////////////////////////////////////////////////////////

MaxEvalFunc = 1000;
NbOuterIter = 25;
GradTOL = 1e-6;
XTol    = 1e-6;
Log     = %T;
Plot    = %T;
PlotLabels = %F;
Plot_ElecOnSPh = %F;

//lambda_init = [0.1,0.1];
lambda_init = [];

[x_opt, x_history, ml_history] = optim_dynq(f, df, ineqconstraint, df_ineqconstraint, eqconstraint, df_eqconstraint, x0, lower, upper, NbOuterIter, GradTOL, XTol, MaxEvalFunc, Log);

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
  printf('inequality constraints:'); disp(ineqconstraint(x_opt)');
end
if (eqconstraint~=[]) then
  printf('equality constraints:'); disp(eqconstraint(x_opt)');
end
printf('Objective function value = %f\n', f(x_opt));

if Plot & ~Plot_ElecOnSPh then
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

  // We plot the inequality constraints
  if (ineqconstraint~=[]) then
    Aux = ineqconstraint(x_opt);
    for k=1:length(Aux)
      for i=1:size(x,2)
        for j=1:size(y,2)
          tmp = ineqconstraint([x(i) y(j)]');
          Z_ineq(i,j) = tmp(k);
        end
      end
      xset('fpf',' ');
      contour(x,y,Z_ineq,[0 0]);
    end
  end

  // We plot the equality constraints
  if (eqconstraint~=[]) then
    Aux = eqconstraint(x_opt);
    for k=1:length(Aux)
      for i=1:size(x,2)
        for j=1:size(y,2)
          tmp = eqconstraint([x(i) y(j)]');
          Z_eq(i,j) = tmp(k);
        end
      end
      xset('fpf',' ');
      contour(x,y,Z_eq,[0 0]);
    end
  end

  for i=1:length(x_history)-1
    plot2d(x_history(i)(1), x_history(i)(2), style=-9);
    if PlotLabels then xstring(x_history(i)(1), x_history(i)(2),string(i)); end

    // Plot the move limits
    Rect(1) = ml_history(i)(2)(1) - ml_history(i)(1)(1);
    Rect(2) = ml_history(i)(2)(2) + ml_history(i)(1)(2);
    Rect(3) = 2*ml_history(i)(1)(1);
    Rect(4) = 2*ml_history(i)(1)(2);
    Color   = -3;
    if i==length(x_history) then Color = -11; end
    xrects(Rect,Color);
  end
  plot2d(x_history($)(1), x_history($)(2), style = -9);
  if PlotLabels then xstring(x_history(i)(1), x_history(i)(2),string(length(x_history))); end
  xtitle('DynQ2 demo','x1','x2');
  drawnow;
end

if Plot & Plot_ElecOnSPh then
  X = matrix(x_opt,N,3);
  param3d1(X(:,1),X(:,2),list(X(:,3),-9));
end

if ineqconstraint~=[] then
  scf();
  Y_ineq = [];
  X_ineq = [];
  for i=1:length(x_history)
    Y_ineq(:,i) = ineqconstraint(x_history(i));
    X_ineq(i)   = i;
  end
  for i=1:size(Y_ineq,1)
    plot(X_ineq,Y_ineq(i,:),'k-');
  end
  plot([X_ineq(1) X_ineq($)],[0 0],'r-.');
  xtitle('Evolution of the inequality constraints','It','Overshoot');
end

if eqconstraint~=[] then
  scf();
  Y_eq = [];
  X_eq = [];
  for i=1:length(x_history)
    Y_eq(:,i) = eqconstraint(x_history(i));
    X_eq(i)   = i;
  end
  for i=1:size(Y_eq,1)
    plot(X_eq,Y_eq(i,:),'k-');
  end
  plot([X_eq(1) X_eq($)],[0 0],'r-.');
  xtitle('Evolution of the equality constraints','It','Overshoot');
end

scf();
Y_obj = [];
X_obj = [];
for i=1:length(x_history)
  Y_obj(i) = f(x_history(i));
  X_obj(i) = i;
end
plot(X_obj,Y_obj,'k-');
xtitle('Evolution of the objective function value','It','f');


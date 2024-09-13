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
//deff('y=eqconstraint(x)','y = larg_sm_poly_eq_constr(x)');
//deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');

/////////////////////////////////////////////////////////////////////////

////function [value] = elec_on_sph_obj(x)
////function [diff_bound] = elec_on_sph_constr(x)
//
//N = 10;
//upper = ones(3*N,1);
//lower = -ones(3*N,1);
//x0 = (upper - lower) .* rand(size(upper,1),size(upper,2)) + lower;
//
//deff('y=f(x)','y = elec_on_sph_obj(matrix(x,N,3))');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=eqconstraint(x)','y = elec_on_sph_ineq_constr(matrix(x,N,3))');
//deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');
//deff('y=eqconstraint(x)','y = elec_on_sph_eq_constr(x)');
//deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');

////////////////////////////////////////////////////////////////////////

////function [value] = vv_open_obj(r)
////function [diff_bound] = vv_open_constr(r)
//
//N = 6;
//upper = ones(N,1);
//lower = zeros(N,1);
//x0 = 0.9*ones(N,1);
//
//deff('y=f(x)','y = vv_open_obj(x)');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=ineqconstraint(x)','y = vv_open_ineq_constr(x)');
//deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');
//deff('y=eqconstraint(x)','y = vv_open_eq_constr(x)');
//deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');

///////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////

////function [value] = himmelblau_1_obj(r)
////function [diff_bound] = himmelblau_1_constr(r)
//
//// 1 3 4 4a 5 8 9 10 11 12 13 14 15 16 24
//PbNumber = '13';
//// Pb avec 13!! log
//// Verifier dans le livre:
//// 3
//// 4 -> mettre abs dans la fobj ??
//// 4a pb sous linux et pas sous windows ??
//// 8 ??
//// 9 il existe des contraintes fortes sur les valeurs de xi. Si ces contraintes ne sont pas respectées, le calcul plante !
//// 13 x(6) apparait dans la fobj. A verifier dans le livre
//
//execstr('x0 = himmelblau_'+PbNumber+'_x_init();');
//
//deff('y=f(x)','y = himmelblau_'+PbNumber+'_obj(x)');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=ineqconstraint(x)','y = himmelblau_'+PbNumber+'_ineq_constr(x)');
//deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');
//deff('y=eqconstraint(x)','y = himmelblau'+PbNumber+'_eq_constr(x)');
//deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');

//////////////////////////////////////////////////////////////////////

//function [value] = G1_obj(r)
//function [diff_bound] = G1_constr(r)

// 1 2 3 4 5  6 7 8 9 10 11 12 13
PbNumber = '2';
// Problem:
// 8 -> verifier ce probleme
// 12, 13 -> pb de convergence -> verifier le probleme
execstr('x0 = G'+PbNumber+'_x_init();');

deff('y=f(x)','y = G'+PbNumber+'_obj(x)');
deff('y=df(x)','y = derivative(f,x)''');
deff('y=ineqconstraint(x)','y = G'+PbNumber+'_ineq_constr(x)');
deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');
deff('y=eqconstraint(x)','y = G'+PbNumber+'_eq_constr(x)');
deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');

//////////////////////////////////////////////////////////////////////

//// weldedbeam
//
//x0 = weldedbeam_x_init()';
//
//deff('y=f(x)','y = weldedbeam_obj(x)');
//deff('y=df(x)','y = derivative(f,x)''');
//deff('y=ineqconstraint(x)','y = weldedbeam_ineq_constr(x)');
//deff('y=df_ineqconstraint(x)','y = derivative(ineqconstraint,x)''');
//deff('y=eqconstraint(x)','y = weldedbeam_eq_constr(x)');
//deff('y=df_eqconstraint(x)','y = derivative(eqconstraint,x)''');

/////////////////////////////////////////////////////////////////////

ItMX    = 500000;
NbLoop  = 10;
Log     = %T;
DispNum = %T;
r_p     = 0.5; // 0.05
//eta     = 2;
eta     = [];
r_max   = 10000;
//r_p     = []; // 0.05
//eta     = [];
//r_max   = [];

Plot = %F;
Plot_SMPoly = %F;
Plot_ElecOnSPh = %F;

g_defined = isdef('ineqconstraint');
h_defined = isdef('eqconstraint');

[x_opt, x_history] = optim_sumt(f, df, ...
                                ineqconstraint, df_ineqconstraint, ...
                                eqconstraint, df_eqconstraint, x0, ItMX, NbLoop, Log, r_p, eta, r_max);

printf('Solution:\n');
//printf('optimal parameters:'); disp(x_opt')
printf('objective function value = %f\n', f(x_opt)); 
if g_defined then
  printf('inequality constraint values :'); disp(max(ineqconstraint(x_opt))');
end
if h_defined then
  printf('equality constraint values :'); disp(max(eqconstraint(x_opt))');
end

printf('\n----------\n\n');

printf('Solution for the initial solution:\n');
//printf('optimal parameters:'); disp(x0')
printf('objective function value = %f\n', f(x0)); 
if g_defined then
  printf('inequality constraint values :'); disp(max(ineqconstraint(x0))');
end
if h_defined then
  printf('equality constraint values :'); disp(max(eqconstraint(x0))');
end

printf('\n----------\n\n');

printf('Solution for each penalty coefficient:\n');
for i=1:length(x_history)
  //printf('Step %d: optimal parameters:', i); disp(x_history(i)')
  printf('Step %d: objective function value = %f\n', i, f(x_history(i))); 
  if g_defined then
    printf('Step %d: inequality constraint values :', i); disp(max(ineqconstraint(x_history(i)))');  
  end
  if h_defined then
    printf('Step %d: equality constraint values :', i); disp(max(eqconstraint(x_history(i)))');  
  end
end

if Plot then
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
end

if Plot_SMPoly then
  X_Plot = x_opt(1:N) .* cos(x_opt(N+1:$));
  Y_Plot = x_opt(1:N) .* sin(x_opt(N+1:$));
  scf();
  plot(X_Plot,Y_Plot,'k-');
end

if Plot_ElecOnSPh then
  X = matrix(x_opt,N,3);
  param3d1(X(:,1),X(:,2),list(X(:,3),-9));
end


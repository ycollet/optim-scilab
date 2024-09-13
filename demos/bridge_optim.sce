//getd('Functions/FEM/femtruss');
global EvalFObj;
global Weight;

UseDiffCode = %F;
UseSciad    = %F;

Weight = 0;

x0 = [1.0; 1.0; 2.0; 1.0; 3.0; 1.0];

function [t,p,e,A,E,rho,F] = bridge_optim(x)
// t : donnees en lignes regroupees par 2 (connection extremite 1 - connection extremite 2)
// p : autant de lignes que de noeuds. Sur chaque ligne n on retrouve les coordonnees 2d du noeud n
// e : liste des noeuds d'appui
// A : liste des sections des elements. Vecteur contenant autant d'entrees que d'elements
// E : liste des modules d'elasticite des elements. Vecteur contenant autant d'entrees que d'elements
// rho : liste des masses volumiques des elements. Vecteur contenant autant d'entrees que d'elements
// F   : liste des forces appliquees aux noeuds. Vecteur colonne comprenant la coordonnes X du noeud 1 en premiere ligne, la coordonnees Y du noeud 1 en seconde 
//       ligne, etc ...

t = [1, 2; // Connectivite element 1
     2, 3; // Connectivite element 2
     3, 4; // Connectivite element 3
     4, 5; // Connectivite element 4
     6, 7; // Connectivite element 5
     7, 8; // Connectivite element 6
     1, 6; // Connectivite element 7
     8, 5; // Connectivite element 8
     6, 2; // Connectivite element 9
     7, 3; // Connectivite element 10
     8, 4; // Connectivite element 11
     6, 3; // Connectivite element 12
     7, 2; // Connectivite element 13
     7, 4; // Connectivite element 14
     8, 3  // Connectivite element 15
    ];

//      6    7    8
//      +----+----+
//     /|\  /|\  /|\  
//    / | /\ | /\ | \ 
//   /  |/  \|/  \|  \
// -+---+----+----+---+-
//  |1  2    3    4   |5

p = [0.0, 0.0; // Noeud 1
     1.0, 0.0; // Noeud 2
     2.0, 0.0; // Noeud 3
     3.0, 0.0; // Noeud 4
     4.0, 0.0; // Noeud 5
     x(1), x(2); // Noeud 6
     x(3), x(4); // Noeud 7
     x(5), x(6)  // Noeud 8
    ];


e = [
     [1, 1] .* localise(5) // On localise les positions du noeud 6 dans les matrices et on immobilise les deux degrees de liberte du noeud 5
     [1, 1] .* localise(1) // On localise les positions du noeud 6 dans les matrices et on immobilise les deux degrees de liberte du noeud 2
    ]; 

A   = ones(1,size(t,1)) * 25e-4; // sections des elements
E   = ones(1,size(t,1)) * 210e9; // module d'elasticite des elements 
rho = ones(1,size(t,1)) * 7.8e3; // masse volumique 

F = [0;  0     // noeud 1
     0; -3.4e5 // noeud 2
     0; -3.4e5 // noeud 3
     0; -3.4e5 // noeud 4
     0;  0     // noeud 5
     0;  0     // noeud 6
     0;  0     // noeud 7
     0;  0     // noeud 8
    ];
endfunction

function y = fobj_truss(x)
global EvalFObj;
global Weight;
k1 = 1.0;
k2 = 0.005;
EvalFObj = EvalFObj + 1;
[t,p,e,A,E,rho,F] = bridge_optim(x);
deff('[t,p,e,A,E,rho,F] = internal_bridge_optim()','[t,p,e,A,E,rho,F] = bridge_optim(x)');
[U,Force,R,T_period,Phi]= femtruss(internal_bridge_optim, %F, %F);
// First objective: minimize the deformation at nodes 2, 3, 4
n = length(U);
U_aux = matrix([U(1:2:n),U(2:2:n)],size(p,1),size(p,2));
Deformation = sqrt(sum(U_aux(2,:).^2 + U_aux(3,:).^2 + U_aux(4,:).^2));
// Second objective: minimize the total length of the bars
Cost = 0;
for i=1:size(t,1)
  Cost = Cost + sum((p(t(i,1),:) - p(t(i,2),:)).^2);
end
y = (1-Weight)*Deformation*k1 + Weight*Cost*k2;
endfunction

function dy = dfobj_truss(x)
ind = 0;
if UseDiffCode then
  printf('Computing directions\n');
  printf('DEBUG: taille vecteur x = %d\n', length(x));
  direction = eye(length(x),length(x));
  printf('Computing the derivatives\n');
  for i=1:length(x)
    x_point = der(x(:),direction(:,i));
    printf('direction %d / %d\n', i, length(x));
    res = fobj_truss(x_point);
    dy(i) = res.v;
  end
elseif UseSciad then
  printf('Computing the graph - part 1\n');
  x_ad = ad_ivar(x);
  printf('Computing the graph - part 2\n');
  cg = fobj_truss(x_ad);
  printf('Computing the objective function value\n');
  f = ad_value(cg);
  direction = eye(length(x),length(x));
  for i=1:length(x)
    printf('Computing direction %d / %d\n', i, length(x));
    //dy(i,1) = ad_jacobian_num(cg,length(x),direction(:,i)); // ad_gradient()
    dy(i,1) = ad_gradient(cg,length(x)); // ad_gradient()
  end
else
  dy = derivative(fobj_truss,x)';
end
endfunction

function y = constr_truss(x)
y(1,1) = x(2) - 2;
y(2,1) = x(4) - 2;
y(3,1) = x(6) - 2;
endfunction

function dy = dconstr_truss(x)
dy = derivative(constr_truss,x)';
endfunction

function [y, dy, ind] = optim_fobj_truss(x, ind)
y  = fobj_truss(x);
dy = derivative(fobj_truss,x);
endfunction

function plot_fobj_truss(x)
[t,p,e,A,E,rho,F] = bridge_optim(x);
deff('[t,p,e,A,E,rho,F] = internal_bridge_optim()','[t,p,e,A,E,rho,F] = bridge_optim(x)');
[U,Force,R,T_period,Phi]= femtruss(internal_bridge_optim, %F, %F);
plotdeforme(U,p,t,10);
endfunction

// Parameters for LFOP
delta_t    = 0.5;  // the temporal step size for the resolution of the dynamical system
delta_step = 0.1;  // maximal allowable step size
GradTOL    = 1e-5; // we stop the algorithm if the gradient is below the GradTOL level.
delta_inc  = 0.01; // increase factor of the dilatation coeff of delta_t
m          = 3;    // id ak+1'.ak is non positive during m iterations then delta_t = delta_t / 2
                     // and we restart from (xk+xk+1) / 2
p_start    = 1.1;  // starting value of the dilatation coefficient: 1.01 by default
MaxEvalFunc = 400;
Algorithm = 'qn'; // 'qn', 'gc', 'nd' -> Ne marche qu'avec 'qn' (quasi-newton). Pour les autres, on obtient rapidement une structure mal
                  // conditionnée
Log = %F;
// Conjugate gradient
h            = 0.001;  // initial step length in line search
ItMX         = 100;    // maximum number of descent steps
ls_ItMX      = 100;    // maximum number of line search steps
lsm  = ls_dicho; // ls_newton = newton line search
                 // ls_secant = secant line search
                 // ls_goldsect = golden section line search
                 // ls_dicho = dichotomy line search
                 // ls_backtrack = backtracking line search
                 // ls_polynom = polynom identification line search
cgm          = 'pr';   // fr = fletcher-reeves
                       // pr = polak-ribiere
                       // hs = hestenes-stiefel
                       // pw = powel
                       // dy = Dai-Yuan
TOL          = 1.0e-6; // accuracy for convergence test (minimum)
StepTOL      = 1.0e-6; // accuracy for convergence test - size of the step 
XTOL         = 1.0e-6; // accuracy for convergence test - improvement on x between two iterations
SubTOL       = 1.0e-6; // accuracy for convergence test (line search)
Restart      = 10;     // restart cg every Restart iterations
NbPtsToCompute  = 9;
MaxWeight       = 0.5;
delta_ml        = 0.1;
MaxMinStepCount = 10;

EvalFObj = 0;
Solutions = list();
Solutions_Cost = [];
Solutions_Deformation = [];
// For ETOP
delta_t_etop = 0.5;
TOL_etop     = 1.0e-12; // accuracy for convergence test (minimum)
XTOL_etop    = 1.0e-12; // accuracy for convergence test - improvement on x between two iterations
ItMX_etop    = 100;
deltax_max   = 0.05;

printf('initial solution:'); disp(x0');
printf('initial objective function value = %f\n',fobj_truss(x0));

for i=0:NbPtsToCompute-1
  // Weigh will vary from 0 to MaxWeight
  printf('Optimisation of point %d / %d\n', i+1, NbPtsToCompute+1);
  Weight = MaxWeight*i/(NbPtsToCompute-1);
  //[x_opt, x_history] = optim_LFOP(dfobj_truss, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, %F, p_start);
  //[x_opt, x_history] = optim_cg(x0, fobj_truss, dfobj_truss, h, Log, ItMX, ls_ItMX, lsm, cgm, TOL, StepTOL, XTOL, SubTOL, Restart);
  //[x_opt, x_history] = optim_bfgs(x0, fobj_truss, dfobj_truss, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart);
  //[x_opt, x_history] = optim_dfp(x0, fobj_truss, dfobj_truss, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart);
  //[x_opt, x_history] = optim_steepest(x0, fobj_truss, dfobj_truss, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL);
  //[f_opt, x_opt] = optim(optim_fobj_truss, x0, algo=Algorithm, iter=MaxEvalFunc);
  [x_opt, x_history] = optim_etop(dfobj_truss, x0, delta_t_etop, ItMX_etop, cgm, deltax_max, XTOL_etop, TOL_etop, Log);

  Solutions($+1) = x_opt;
  
  printf('number of call to objective function: %d\n',EvalFObj);

  printf('Final solution:'); disp(x_opt');
  printf('Final objective function value = %f\n',fobj_truss(x_opt));
  
  Weight = 1; Solutions_Cost($+1)        = fobj_truss(x_opt);
  Weight = 0; Solutions_Deformation($+1) = fobj_truss(x_opt);
end
plot(Solutions_Cost, Solutions_Deformation, 'ro');
xtitle('Pareto front','Length of bars','Deformation');
scf();
plot_fobj_truss(x0);
xtitle('Before optimization','x','y');
scf();
Index = 1;
NbCol = ceil(sqrt(NbPtsToCompute));
for i=1:NbPtsToCompute
  Weight = MaxWeight*i/(NbPtsToCompute-1);
  subplot(NbCol,NbCol,i);
  plot_fobj_truss(Solutions(i));
  xtitle(sprintf('After optimization - Weight = %f', Weight),'x','y');
end

[x_opt, x_history] = optim_slp(fobj_truss, dfobj_truss, constr_truss, dconstr_truss, [], [], x0, ItMX, MaxEvalFunc, delta_ml, ...
                               7*ones(size(x0,1),size(x0,2)), zeros(size(x0,1),size(x0,2)), MaxMinStepCount, Log);

printf('Constr - Final solution:'); disp(x_opt');
printf('Constr - Final objective function value = %f\n',fobj_truss(x_opt));
scf();
drawlater;
plot_fobj_truss(x_opt);
xtitle('after constrained optimization','x','y');
drawnow;

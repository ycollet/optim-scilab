printf('Optimization of a pole structure\n\n');

Log = %F;

global EvalFObj;
global Weight;

// On ne touche pas aux noeuds 1, 2, 16, 17, 18, 19 - les 4 derniers supportents les cables
// les  deux premiers correspondent a la prise du pylone sur la terre

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5

x0 = [2.0; 1.0; ... // OK
      4.0; 1.0; ... // OK
      2.0; 2.0; ... // OK
      4.0; 2.0; ... // OK
      2.0; 3.0; ... // OK
      4.0; 3.0; ... // OK
      2.0; 4.0; ... // OK
      4.0; 4.0; ... // OK
      2.0; 5.0; ... // OK
      4.0; 5.0; ... // OK
      2.0; 6.0; ... // OK
      4.0; 6.0; ... // OK
      3.0; 7.0; ... // OK
      1.0; 5.0; ... // OK
      5.0; 5.0];    // OK

function [t,p,e,A,E,rho,F] = pylone_optim(x)
// t   : donnees en lignes regroupees par 2 (connection extremite 1 - connection extremite 2)
// p   : autant de lignes que de noeuds. Sur chaque ligne n on retrouve les coordonnees 2d du noeud n
// e   : liste des noeuds d'appui
// A   : liste des sections des elements. Vecteur contenant autant d'entrees que d'elements
// E   : liste des modules d'elasticite des elements. Vecteur contenant autant d'entrees que d'elements
// rho : liste des masses volumiques des elements. Vecteur contenant autant d'entrees que d'elements
// F   : liste des forces appliquees aux noeuds. Vecteur colonne comprenant la coordonnes X du noeud 1 en premiere ligne, la coordonnees Y du noeud 1 en seconde 
//       ligne, etc ...

t = [1,  2; // Les barres horizontales
     3,  4;
     5,  6;
     7,  8;
     9,  10;
     11, 12;
     13, 14;
     16, 17;
     17, 9;
     10, 18;
     18, 19;
     20, 11;
     12, 21;
     1,  3; // Les barres verticales
     2,  4;
     3,  5;
     4,  6;
     5,  7;
     6,  8;
     7,  9;
     8,  10;
     9,  11;
     10, 12;
     11, 13;
     12, 14;
     17, 20;
     18, 21;
     1,  4; // Les barres diagonales droites
     3,  6;
     5,  8;
     7,  10;
     9,  12;
     11, 14;
     13, 15;
     17, 11;
     20, 13;
     16, 20;
     10, 21;
     2,  3; // Les barres de diagonales gauches
     4,  5;
     6,  7;
     8,  9;
     10, 11;
     12, 13;
     14, 15;
     19, 21;
     18, 12;
     21, 14;
     9,  20
     ];

p = [2.0,   0.0;   // Noeud 1
     4.0,   0.0;   // Noeud 2
     x(1),  x(2);  // Noeud 3
     x(3),  x(4);  // Noeud 4
     x(5),  x(6);  // Noeud 5
     x(7),  x(8);  // Noeud 6
     x(9),  x(10); // Noeud 7
     x(11), x(12); // Noeud 8
     x(13), x(14); // Noeud 9
     x(15), x(16); // Noeud 10
     x(17), x(18); // Noeud 11
     x(19), x(20); // Noeud 12
     x(21), x(22); // Noeud 13
     x(23), x(24); // Noeud 14
     x(25), x(26); // Noeud 15
     0.0,   4.0;   // Noeud 16
     1.0,   4.0;   // Noeud 17
     5.0,   4.0;   // Noeud 18
     6.0,   4.0;   // Noeud 19
     x(27), x(28); // Noeud 20
     x(29), x(30)  // Noeud 21
    ];

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5


e = [
     [1, 1] .* localise(1) // On localise les positions du noeud 1 dans les matrices et on immobilise les deux degrees de liberte du noeud 1 et 2
     [1, 1] .* localise(2) // On localise les positions du noeud 2 dans les matrices et on immobilise les deux degrees de liberte du noeud 1 et 2
    ]; 

A   = ones(1,size(t,1)) * 25e-4; // sections des elements
E   = ones(1,size(t,1)) * 210e9; // module d'elasticite des elements 
rho = ones(1,size(t,1)) * 7.8e3; // masse volumique 

// Le pylone soutient 4 cables de 100 Kg chacun
F = [0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; -1e5;
     0.0; -1e5;
     0.0; -1e5;
     0.0; -1e5;
     0.0; 0.0;
     0.0; 0.0
    ];
endfunction

function y = fobj_truss(x)
global EvalFObj;
global Weight;
k1 = 1.0;
k2 = 0.005;
EvalFObj = EvalFObj + 1;
[t,p,e,A,E,rho,F] = pylone_optim(x);
deff('[t,p,e,A,E,rho,F] = internal_pylone_optim()','[t,p,e,A,E,rho,F] = pylone_optim(x)');
[U,Force,R,T_period,Phi]= femtruss(internal_pylone_optim, %F, Log);
// First objective: minimize the deformation at nodes 16, 17, 18, 19
n = length(U);
U_aux = matrix([U(1:2:n),U(2:2:n)],size(p,1),size(p,2));
Deformation = sqrt(sum(U_aux(16,:).^2 + U_aux(17,:).^2 + U_aux(18,:).^2 + U_aux(19,:).^2));
// Second objective: minimize the total length of the bars
Cost = 0;
for i=1:size(t,1)
  Cost = Cost + sum((p(t(i,1),:) - p(t(i,2),:)).^2);
end
y = (1-Weight)*Deformation*k1 + Weight*Cost*k2;
endfunction

function dy = dfobj_truss(x)
dy = derivative(fobj_truss,x)';
endfunction

function [y, dy, ind] = optim_fobj_truss(x, ind)
y  = fobj_truss(x);
dy = derivative(fobj_truss,x);
endfunction

function plot_fobj_truss(x)
[t,p,e,A,E,rho,F] = pylone_optim(x);
deff('[t,p,e,A,E,rho,F] = internal_pylone_optim()','[t,p,e,A,E,rho,F] = pylone_optim(x)');
[U,Force,R,T_period,Phi]= femtruss(internal_pylone_optim, %F, Log);
plotdeforme(U,p,t,10);
endfunction

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5

function y=constr_truss(x)
y(1,1) = 4.0   - x(28);
y(2,1) = 4.0   - x(30);
y(3,1) = x(27) - x(17);
y(4,1) = 1.0   - x(13);
y(5,1) = x(19) - x(29);
y(6,1) = x(15) - 5.0;

// Constraints related to the length of the bars
// They must be longer than 0.5
[t,p,e,A,E,rho,F] = pylone_optim(x);
Index = 7;
for i=1:size(t,1)
  if ((~isempty(find([1,2,16,17,18,19]==t(i,1))))&(~isempty(find([1,2,16,17,18,19]==t(i,2))))) then
    // If a constraint involve two fixed points, we skip it. Because,otherwise, the gradient of the constraint will be null.
    continue;
  end
  y(Index,1) = - sqrt( (p(t(i,1),1) - p(t(i,2),1))^2 + (p(t(i,1),2) - p(t(i,2),2))^2) + 0.5;
  Index = Index + 1;
end

// Constraints related to the length of the bars
// They must be smaller than 3
for i=1:size(t,1)
  if ((~isempty(find([1,2,16,17,18,19]==t(i,1))))&(~isempty(find([1,2,16,17,18,19]==t(i,2))))) then
    // If a constraint involve two fixed points, we skip it. Because,otherwise, the gradient of the constraint will be null.
    continue;
  end
  y(Index,1) = sqrt( (p(t(i,1),1) - p(t(i,2),1))^2 + (p(t(i,1),2) - p(t(i,2),2))^2) - 3;
  Index = Index + 1;
end
endfunction

function dy=dconstr_truss(x)
dy = derivative(constr_truss,x)';
endfunction

/////////////////////////////////
// Symetric objective function //
/////////////////////////////////

x0_sym = [2.0; 1.0; ...
          2.0; 2.0; ...
          2.0; 3.0; ...
          2.0; 4.0; ...
          2.0; 5.0; ...
          2.0; 6.0; ...
          2.5; 7.0; ...
          1.0; 5.0];

function [t,p,e,A,E,rho,F] = pylone_optim_sym(x)
// t   : donnees en lignes regroupees par 2 (connection extremite 1 - connection extremite 2)
// p   : autant de lignes que de noeuds. Sur chaque ligne n on retrouve les coordonnees 2d du noeud n
// e   : liste des noeuds d'appui
// A   : liste des sections des elements. Vecteur contenant autant d'entrees que d'elements
// E   : liste des modules d'elasticite des elements. Vecteur contenant autant d'entrees que d'elements
// rho : liste des masses volumiques des elements. Vecteur contenant autant d'entrees que d'elements
// F   : liste des forces appliquees aux noeuds. Vecteur colonne comprenant la coordonnes X du noeud 1 en premiere ligne, la coordonnees Y du noeud 1 en seconde 
//       ligne, etc ...

t = [1,  2; // Les barres horizontales
     3,  4;
     5,  6;
     7,  8;
     9,  10;
     11, 12;
     13, 14;
     16, 17;
     17, 9;
     10, 18;
     18, 19;
     20, 11;
     12, 21;
     1,  3; // Les barres verticales
     2,  4;
     3,  5;
     4,  6;
     5,  7;
     6,  8;
     7,  9;
     8,  10;
     9,  11;
     10, 12;
     11, 13;
     12, 14;
     17, 20;
     18, 21;
     1,  4; // Les barres diagonales droites
     3,  6;
     5,  8;
     7,  10;
     9,  12;
     11, 14;
     13, 15;
     17, 11;
     20, 13;
     16, 20;
     10, 21;
     2,  3; // Les barres de diagonales gauches
     4,  5;
     6,  7;
     8,  9;
     10, 11;
     12, 13;
     14, 15;
     19, 21;
     18, 12;
     21, 14;
     9,  20
     ];

p = [2.0,       0.0;   // Noeud 1
     3.0,       0.0;   // Noeud 2
     x(1),      x(2);  // Noeud 3
     5 - x(1),  x(2);  // Noeud 4
     x(3),      x(4);  // Noeud 5
     5 - x(3),  x(4);  // Noeud 6
     x(5),      x(6);  // Noeud 7
     5 - x(5),  x(6);  // Noeud 8
     x(7),      x(8);  // Noeud 9
     5 - x(7),  x(8);  // Noeud 10
     x(9),      x(10); // Noeud 11
     5 - x(9),  x(10); // Noeud 12
     x(11),     x(12); // Noeud 13
     5 - x(11), x(12); // Noeud 14
     x(13),     x(14); // Noeud 15
     0.0,       4.0;   // Noeud 16
     1.0,       4.0;   // Noeud 17
     4.0,       4.0;   // Noeud 18
     5.0,       4.0;   // Noeud 19
     x(15),     x(16); // Noeud 20
     5 - x(15), x(16)  // Noeud 21
    ];

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5


e = [
     [1, 1] .* localise(1) // On localise les positions du noeud 1 dans les matrices et on immobilise les deux degrees de liberte du noeud 1 et 2
     [1, 1] .* localise(2) // On localise les positions du noeud 2 dans les matrices et on immobilise les deux degrees de liberte du noeud 1 et 2
    ]; 

A   = ones(1,size(t,1)) * 25e-4; // sections des elements
E   = ones(1,size(t,1)) * 210e9; // module d'elasticite des elements 
rho = ones(1,size(t,1)) * 7.8e3; // masse volumique 

// Le pylone soutient 4 cables de 100 Kg chacun
F = [0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; 0.0;
     0.0; -1e5;
     0.0; -1e5;
     0.0; -1e5;
     0.0; -1e5;
     0.0; 0.0;
     0.0; 0.0
    ];
endfunction

function y = fobj_truss_sym(x)
global EvalFObj;
global Weight;
k1 = 1.0;
k2 = 0.005;
EvalFObj = EvalFObj + 1;
[t,p,e,A,E,rho,F] = pylone_optim_sym(x);
deff('[t,p,e,A,E,rho,F] = internal_pylone_optim_sym()','[t,p,e,A,E,rho,F] = pylone_optim_sym(x)');
[U,Force,R,T_period,Phi]= femtruss(internal_pylone_optim_sym, %F, Log);
// First objective: minimize the deformation at nodes 16, 17, 18, 19
n = length(U);
U_aux = matrix([U(1:2:n),U(2:2:n)],size(p,1),size(p,2));
Deformation = sqrt(sum(U_aux(16,:).^2 + U_aux(17,:).^2 + U_aux(18,:).^2 + U_aux(19,:).^2));
// Second objective: minimize the total length of the bars
Cost = 0;
for i=1:size(t,1)
  Cost = Cost + sum((p(t(i,1),:) - p(t(i,2),:)).^2);
end
y = (1-Weight)*Deformation*k1 + Weight*Cost*k2;
endfunction

function dy = dfobj_truss_sym(x)
dy = derivative(fobj_truss_sym,x)';
endfunction

function [y, dy, ind] = optim_fobj_truss_sym(x, ind)
y  = fobj_truss_sym(x);
dy = derivative(fobj_truss_sym,x);
endfunction

function plot_fobj_truss_sym(x)
[t,p,e,A,E,rho,F] = pylone_optim_sym(x);
deff('[t,p,e,A,E,rho,F] = internal_pylone_optim_sym()','[t,p,e,A,E,rho,F] = pylone_optim_sym(x)');
[U,Force,R,T_period,Phi]= femtruss(internal_pylone_optim_sym, %F, Log);
plotdeforme(U,p,t,10);
endfunction

// 7                    + (15)
// 6               + (13)    + (14)
// 5        + (20) + (11)    + (12) + (21)
// 4 + (16) + (17) + (9)     + (10) + (18) + (19)
// 3               + (7)     + (8)
// 2               + (5)     + (6)
// 1               + (3)     + (4)
// 0               + (1)     + (2)
//   0      1      2         3      4      5


//p = [2.0,       0.0;   // 2.0,   0.0;   // Noeud 1
//     4.0,       0.0;   // 4.0,   0.0;   // Noeud 2
//     x(1),      x(2);  // x(1),  x(2);  // Noeud 3
//     5 - x(1),  x(2);  // x(3),  x(4);  // Noeud 4
//     x(3),      x(4);  // x(5),  x(6);  // Noeud 5
//     5 - x(3),  x(4);  // x(7),  x(8);  // Noeud 6
//     x(5),      x(6);  // x(9),  x(10); // Noeud 7
//     5 - x(5),  x(6);  // x(11), x(12); // Noeud 8
//     x(7),      x(8);  // x(13), x(14); // Noeud 9
//     5 - x(7),  x(8);  // x(15), x(16); // Noeud 10
//     x(9),      x(10); // x(17), x(18); // Noeud 11
//     5 - x(9),  x(10); // x(19), x(20); // Noeud 12
//     x(11),     x(12); // x(21), x(22); // Noeud 13
//     5 - x(11), x(12); // x(23), x(24); // Noeud 14
//     x(13),     x(14); // x(25), x(26); // Noeud 15
//     0.0,       4.0;   // 0.0,   4.0;   // Noeud 16
//     1.0,       4.0;   // 1.0,   4.0;   // Noeud 17
//     5.0,       4.0;   // 5.0,   4.0;   // Noeud 18
//     6.0,       4.0;   // 6.0,   4.0;   // Noeud 19
//     x(15),     x(16); // x(27), x(28); // Noeud 20
//     5 - x(15), x(16)  // x(29), x(30)  // Noeud 21
//    ];

function y=constr_truss_sym(x)
y(1,1) = 4.0   - x(16);
y(2,1) = x(15) - x(9);
y(3,1) = 1.0   - x(13);
y(4,1) = (5 - x(9)) - (5 - x(15));
y(5,1) = 5 - x(7) - 5.0;

// Constraints related to the length of the bars
// They must be longer than 0.5
[t,p,e,A,E,rho,F] = pylone_optim_sym(x);
Index = 6;
for i=1:size(t,1)
  if ((~isempty(find([1,2,16,17,18,19]==t(i,1))))&(~isempty(find([1,2,16,17,18,19]==t(i,2))))) then
    // If a constraint involve two fixed points, we skip it. Because,otherwise, the gradient of the constraint will be null.
    continue;
  end
  y(Index,1) = - sqrt( (p(t(i,1),1) - p(t(i,2),1))^2 + (p(t(i,1),2) - p(t(i,2),2))^2) + 0.5;
  Index = Index + 1;
end
endfunction

function dy=dconstr_truss_sym(x)
dy = derivative(constr_truss_sym,x)';
endfunction

// Parameters for LFOP and LFOPc
delta_t    = 0.05;  // the temporal step size for the resolution of the dynamical system YC:0.5 avant
delta_step = 0.1;  // maximal allowable step size
GradTOL    = 1e-5; // we stop the algorithm if the gradient is below the GradTOL level.
delta_inc  = 0.01; // increase factor of the dilatation coeff of delta_t
m          = 3;    // id ak+1'.ak is non positive during m iterations then delta_t = delta_t / 2
                   // and we restart from (xk+xk+1) / 2
p_start    = 1.1;  // starting value of the dilatation coefficient: 1.01 by default
MaxEvalFunc = 1000;
Algorithm = 'qn'; // 'qn', 'gc', 'nd' -> Ne marche qu'avec 'qn' (quasi-newton). Pour les autres, on obtient rapidement une structure mal
                  // conditionnée
// Conjugate gradient
h            = 0.001; // initial step length in line search
ItMX         = 100;   // maximum number of descent steps
ls_ItMX      = 100;   // maximum number of line search steps
lsm  = ls_dicho; // ls_newton    = newton line search
                 // ls_secant    = secant line search
                 // ls_goldsect  = golden section line search
                 // ls_dicho     = dichotomy line search
                 // ls_backtrack = backtracking line search
                 // ls_polynom   = polynom identification line search
cgm = 'fr'; // fr = fletcher-reeves
            // pr = polak-ribiere
            // hs = hestenes-stiefel
            // pw = powel
            // dy = Dai-Yuan
TOL            = 1.0e-6; // accuracy for convergence test (minimum)
StepTOL        = 1.0e-6; // accuracy for convergence test - size of the step 
XTOL           = 1.0e-6; // accuracy for convergence test - improvement on x between two iterations
SubTOL         = 1.0e-6; // accuracy for convergence test (line search)
Restart        = 10;     // restart cg every Restart iterations
NbPtsToCompute = 1;
MaxWeight      = 0.5;
Constraints    = %T;
MaxMinStepCount = 3;
delta_ml        = 0.2; // 0.2 avant
NbOuterIter     = 5;

UseSymetric = %T; // Use the symetric formulation of the problem
rand_ampl   = 0.001; // Amplitude of the value value we add to x0 -> to try ta have several solutions
// SUMT
NbLoop = 3;
r_p    = 2.0;
eta    = 2.0;
r_max  = 2.0;
// Feasdir
SubTOL = 1e-2;
Theta_0 = 1.0;
epsilon = -0.01;
// SQSD
rho = 0.1;
// etop
delta_t_etop = 0.5;
TOL_etop     = 1.0e-12; // accuracy for convergence test (minimum)
XTOL_etop    = 1.0e-12; // accuracy for convergence test - improvement on x between two iterations
ItMX_etop    = 100;
deltax_max   = 0.05;

EvalFObj = 0;
Solutions = list();
Solutions_Cost = [];
Solutions_Deformation = [];

if (UseSymetric) then
  x0_sym = x0_sym + rand_ampl*rand(size(x0_sym,1),size(x0_sym,2));
else
  x0 = x0 + rand_ampl*rand(size(x0,1),size(x0,2));
end

if (~Constraints) then
  printf('optimization without constraints\n');
  if (UseSymetric) then
    printf('using symetric parametrization\n');
    Weight = 1; Solutions_Cost($+1)        = fobj_truss_sym(x0_sym);
    Weight = 0; Solutions_Deformation($+1) = fobj_truss_sym(x0_sym);
    printf('initial solution:'); disp(x0_sym');
  else
    printf('using non-symetric parametrization\n');
    Weight = 1; Solutions_Cost($+1)        = fobj_truss(x0);
    Weight = 0; Solutions_Deformation($+1) = fobj_truss(x0);
    printf('initial solution:'); disp(x0');
  end

  printf('initial objective function value = %f\n',Solutions_Deformation($));

  for i=0:NbPtsToCompute-1
    // Weigh will vary from 0 to MaxWeight
    printf('Optimisation of point %d / %d\n', i+1, NbPtsToCompute+1);
    if (NbPtsToCompute==1) then
      Weight = 0;
    else
      Weight = MaxWeight*i/(NbPtsToCompute-1);
    end
    
    if (UseSymetric) then
      [x_opt, x_history] = optim_etop(dfobj_truss_sym, x0_sym, delta_t_etop, ItMX_etop, cgm, deltax_max, XTOL_etop, TOL_etop, Log);
      //[x_opt, x_history] = optim_sqsd(fobj_truss_sym, dfobj_truss_sym, x0_sym, ItMX, XTOL, TOL, rho, Log);
      //[x_opt, x_history] = optim_LFOP(dfobj_truss_sym, x0_sym, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
      //[x_opt, x_history] = optim_cg(x0_sym, fobj_truss_sym, dfobj_truss_sym, h, Log, ItMX, ls_ItMX, lsm, cgm, TOL, StepTOL, XTOL, SubTOL, Restart);
      //[x_opt, x_history] = optim_bfgs(x0_sym, fobj_truss_sym, dfobj_truss_sym, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart);
      //[x_opt, x_history] = optim_dfp(x0_sym, fobj_truss_sym, dfobj_truss_sym, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart);
      //[x_opt, x_history] = optim_steepest(x0_sym, fobj_truss_sym, dfobj_truss_sym, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL);
      //[f_opt, x_opt] = optim(optim_fobj_truss_sym, x0_sym, algo=Algorithm, iter=MaxEvalFunc);

//      // The method optim_etop must be used after another optimization method
//      printf('ETOP: Final objective function value = %f\n',fobj_truss_sym(x_opt));
//      printf('number of call to objective function: %d\n',EvalFObj);
//      Weight_old = Weight;
//      Weight = 1; Solutions_Cost($+1)        = fobj_truss_sym(x_opt);
//      Weight = 0; Solutions_Deformation($+1) = fobj_truss_sym(x_opt);
//      Weight = Weight_old;
//      scf();
//      plot_fobj_truss_sym(x_opt);
//      xtitle('Before call to etop','x','y');
//      [x_opt, x_history] = optim_etop(dfobj_truss_sym, x_opt, delta_t_etop, ItMX_etop, cgm, XTOL_etop, TOL_etop, Log);
//      plot_fobj_truss_sym(x_opt);
    else
      [x_opt, x_history] = optim_etop(dfobj_truss, x0, delta_t_etop, ItMX_etop, cgm, XTOL_etop, TOL_etop, Log);
      //[x_opt, x_history] = optim_sqsd(fobj_truss, dfobj_truss, x0, ItMX, XTOL, TOL, rho, Log);
      //[x_opt, x_history] = optim_LFOP(dfobj_truss, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
      //[x_opt, x_history] = optim_cg(x0, fobj_truss, dfobj_truss, h, Log, ItMX, ls_ItMX, lsm, cgm, TOL, StepTOL, XTOL, SubTOL, Restart);
      //[x_opt, x_history] = optim_bfgs(x0, fobj_truss, dfobj_truss, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart);
      //[x_opt, x_history] = optim_dfp(x0, fobj_truss, dfobj_truss, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart);
      //[x_opt, x_history] = optim_steepest(x0, fobj_truss, dfobj_truss, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL);
      //[f_opt, x_opt] = optim(optim_fobj_truss, x0, algo=Algorithm, iter=MaxEvalFunc);
  
//      // The method optim_etop must be used after another optimization method
//      printf('ETOP: Final objective function value = %f\n',fobj_truss(x_opt));
//      printf('number of call to objective function: %d\n',EvalFObj);
//      Weight_old = Weight;
//      Weight = 1; Solutions_Cost($+1)        = fobj_truss(x_opt);
//      Weight = 0; Solutions_Deformation($+1) = fobj_truss(x_opt);
//      Weight = Weight_old;
//      scf();
//      plot_fobj_truss(x_opt);
//      xtitle('Before call to etop','x','y');
//      [x_opt, x_history] = optim_etop(dfobj_truss, x_opt, delta_t_etop, ItMX_etop, cgm, XTOL_etop, TOL_etop, Log);
//      plot_fobj_truss(x_opt);
    end

    Solutions($+1) = x_opt;

    printf('Final solution:'); disp(x_opt');
    
    if (UseSymetric) then
      printf('Final objective function value = %f\n',fobj_truss_sym(x_opt));

      Weight = 1; Solutions_Cost($+1)        = fobj_truss_sym(x_opt);
      Weight = 0; Solutions_Deformation($+1) = fobj_truss_sym(x_opt);
    else
      printf('Final objective function value = %f\n',fobj_truss(x_opt));

      Weight = 1; Solutions_Cost($+1)        = fobj_truss(x_opt);
      Weight = 0; Solutions_Deformation($+1) = fobj_truss(x_opt);
    end
    printf('number of call to objective function: %d\n',EvalFObj);
  end

  scf();
  plot(Solutions_Cost, Solutions_Deformation, 'ro');
  xtitle('Pareto front','Length of bars','Deformation');
  
  scf();
  if (UseSymetric) then
    plot_fobj_truss_sym(x0_sym);
  else
    plot_fobj_truss(x0);
  end
  xtitle('Before optimization','x','y');
  scf();
  if (NbPtsToCompute==1) then
    if (UseSymetric) then
      plot_fobj_truss_sym(Solutions(1));
    else
      plot_fobj_truss(Solutions(1));
    end
    xtitle(sprintf('After optimization - Weight = %f', Weight),'x','y');
  else
    Index = 1;
    NbCol = ceil(sqrt(NbPtsToCompute));
    for i=1:NbPtsToCompute
      Weight = MaxWeight*i/(NbPtsToCompute-1);
      subplot(NbCol,NbCol,i);
      if (UseSymetric) then
        plot_fobj_truss_sym(Solutions(i));
      else
        plot_fobj_truss(Solutions(i));
      end
      xtitle(sprintf('After optimization - Weight = %f', Weight),'x','y');
    end
  end
  X_test = [];
  Index  = 1;
  IndexList = [];
  for i=1:length(x_history)
    if (typeof(x_history(i))=='list') then
      IndexList($+1) = Index;
      for j=1:length(x_history(i))
        if (UseSymetric) then
          X_test($+1) = fobj_truss_sym(x_history(i)(j));
        else
          X_test($+1) = fobj_truss(x_history(i)(j));
        end
        Index = Index + 1;
      end
    else
      if (UseSymetric) then
        X_test($+1) = fobj_truss_sym(x_history(i));
      else
        X_test($+1) = fobj_truss(x_history(i));
      end
    end
  end
  scf();
  plot(X_test);
  if (~isempty(IndexList)) then
    plot(IndexList, X_test(IndexList),'ro');
  end
  xtitle('Evolution of the objective function','Iteration','FObj');
else
  printf('optimization without constraints\n');
  if (UseSymetric) then
    printf('using symetric parametrization\n');
    printf('initial solution:'); disp(x0_sym');
    printf('initial objective function value = %f\n',fobj_truss_sym(x0_sym));
    printf('initial constraints value = '); disp(constr_truss_sym(x0_sym));
  else
    printf('using non-symetric parametrization\n');
    printf('initial solution:'); disp(x0');
    printf('initial objective function value = %f\n',fobj_truss(x0));
    printf('initial constraints value = '); disp(constr_truss(x0));
  end
  
  if (UseSymetric) then
    Weight = 1; Solutions_Cost($+1)        = fobj_truss_sym(x0_sym);
    Weight = 0; Solutions_Deformation($+1) = fobj_truss_sym(x0_sym);
    printf('initial solution:'); disp(x0_sym');
  else
    Weight = 1; Solutions_Cost($+1)        = fobj_truss(x0);
    Weight = 0; Solutions_Deformation($+1) = fobj_truss(x0);
    printf('initial solution:'); disp(x0');
  end

  printf('initial objective function value = %f\n',Solutions_Deformation($));

  for i=0:NbPtsToCompute-1
    // Weigh will vary from 0 to MaxWeight
    printf('Optimisation of point %d / %d\n', i+1, NbPtsToCompute+1);
    if (NbPtsToCompute==1) then
      Weight = 0;
    else
      Weight = MaxWeight*i/(NbPtsToCompute-1);
    end

    if (UseSymetric) then
      upper = 7*ones(length(x0_sym),1);
      lower = zeros(length(x0_sym),1);
      // [x_opt, x_history] = optim_slp(fobj_truss_sym, dfobj_truss_sym, constr_truss_sym, dconstr_truss_sym, [], [], x0_sym, ItMX, MaxEvalFunc, delta_ml, upper, lower, MaxMinStepCount, Log);
      // [x_opt, x_history] = optim_LFOPc(dfobj_truss_sym, constr_truss_sym, dconstr_truss_sym, [], [], x0_sym, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
      // [x_opt, x_history] = optim_sumt(fobj_truss_sym, dfobj_truss_sym, constr_truss_sym, dconstr_truss_sym, [], [], x0_sym, ItMX, NbLoop, Log, r_p, eta, r_max);
      [x_opt, x_history] = optim_dynq(fobj_truss_sym, dfobj_truss_sym, constr_truss_sym, dconstr_truss_sym, [], [], x0_sym, lower, upper, NbOuterIter, GradTOL, 0, MaxEvalFunc, %F);
      // [x_opt, x_history] = optim_mma(fobj_truss_sym, dfobj_truss_sym, constr_truss_sym, dconstr_truss_sym, x0_sym, upper, lower, ItMX, MaxEvalFunc, delta_step, Log);
      // [x_opt, x_history] = optim_feasdir(fobj_truss_sym, dfobj_truss_sym, constr_truss_sym, dconstr_truss_sym, x0_sym, ItMX, MaxEvalFunc, upper, lower, Theta_0, epsilon, h, MaxMinStepCount, SubTOL, ls_ItMX, Log);

      printf('number of call to objective function: %d\n',EvalFObj);

      printf('Final solution: f = %f\n',fobj_truss_sym(x_opt));
      printf('Final constraints value ='); disp(constr_truss_sym(x_opt));

      plot_fobj_truss_sym(x_opt);

      xtitle('Constrained optimization','x1','x2');

      Weight = 1; Solutions_Cost($+1)        = fobj_truss_sym(x_opt);
      Weight = 0; Solutions_Deformation($+1) = fobj_truss_sym(x_opt);
    else
      upper = 7*ones(length(x0),1);
      lower = zeros(length(x0),1);
      // [x_opt, x_history] = optim_slp(fobj_truss, dfobj_truss, constr_truss, dconstr_truss, [], [], x0, ItMX, MaxEvalFunc, delta_ml, upper, lower, MaxMinStepCount, Log);
      // [x_opt, x_history] = optim_LFOPc(dfobj_truss, constr_truss, dconstr_truss, [], [], x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start);
      [x_opt, x_history] = optim_dynq(fobj_truss, dfobj_truss, constr_truss, dconstr_truss, [], [], x0, lower, upper, NbOuterIter, GradTOL, 0, MaxEvalFunc, Log);
      // [x_opt, x_history] = optim_sumt(fobj_truss, dfobj_truss, constr_truss, dconstr_truss, [], [], x0, ItMX, NbLoop, Log, r_p, eta, r_max);
      // [x_opt, x_history] = optim_mma(fobj_truss, dfobj_truss, constr_truss, dconstr_truss, x0, upper, lower, ItMX, MaxEvalFunc, delta_step, Log)
      // [x_opt, x_history] = optim_feasdir(fobj_truss, dfobj_truss, constr_truss, dconstr_truss, x0, ItMX, MaxEvalFunc, upper, lower, Theta_0, epsilon, h, MaxMinStepCount, SubTOL, ls_ItMX, Log);

      printf('number of call to objective function: %d\n',EvalFObj);

      printf('Final solution: f = %f\n',fobj_truss(x_opt));
      printf('Final constraints value ='); disp(constr_truss(x_opt));

      plot_fobj_truss(x_opt);

      xtitle('Constrained optimization','x1','x2');

      Weight = 1; Solutions_Cost($+1)        = fobj_truss(x_opt);
      Weight = 0; Solutions_Deformation($+1) = fobj_truss(x_opt);
    end
  end
  
  X_test = [];
  Index  = 1;
  IndexList = [];
  for i=1:length(x_history)
    if (typeof(x_history(i))=='list') then
      IndexList($+1) = Index;
      for j=1:length(x_history(i))
        if (UseSymetric) then
          X_test($+1) = fobj_truss_sym(x_history(i)(j));
        else
          X_test($+1) = fobj_truss(x_history(i)(j));
        end
        Index = Index + 1;
      end
    else
      if (UseSymetric) then
        X_test($+1) = fobj_truss_sym(x_history(i));
      else
        X_test($+1) = fobj_truss(x_history(i));
      end
    end
  end
  scf();
  if (~isempty(IndexList)) then
    subplot(2,1,1);
  end
  plot(X_test);
  xtitle('All the value of objective functions','Iterations','FObj');
  if (~isempty(IndexList)) then
    subplot(2,1,2);
    plot(IndexList, X_test(IndexList),'ro');
    xtitle('The value of objective functions for each optimization sequences','Iterations','FObj');
  end
end

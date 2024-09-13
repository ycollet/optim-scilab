//funcname = 'paviani'; // Attention: xi~=0 et xi~=10
//funcname = 'dixonprice';
//funcname = 'levy8';// A verifier
//funcname = 'levy2';// A verifier
//funcname = 'holtzmann'; // LBFGS line 77 de lsm
//funcname = 'gen_rosen';
//funcname = 'griewank'; // LBFGS line 77 de lsm
//funcname = 'sphere'; // LBFGS line 77 de lsm
//funcname = 'weierstrass';
//funcname = 'ackley'; // LBFGS line 77 de lsm
//funcname = 'ellipsoid'; // LBFGS  line 77 de lsm
//funcname = 'rotell'; // LBFGS  line 77 de lsm
//funcname = 'abspow'; // LBFGS  line 77 de lsm
//funcname = 'michalewicz'; // LBFGS  line 77 de lsm
//funcname = 'powell'; // LBFGS  line 77 de lsm
//funcname = 'gen_rastrigin'; // LBFGS  line 77 de lsm
//funcname = 'schwefel';
//funcname = 'trid';
funcname = 'zhakarov'; // LBFGS  line 77 de lsm
//funcname = 'zhufu'; // pb d'estimation de derive par derivative

// Conjugate gradient method
cgm = 'pr';   // fr = fletcher-reeves
              // pr = polak-ribiere
              // hs = hestenes-stiefel
              // pw = powel
              // dy = Dai-Yuan

// General parameters for optimization methods
ItMX     = 500;         // maximum number of descent steps
ls_ItMX  = 100;        // maximum number of line search iteration
Restart  = 10;         // restart cg every Restart iterations
TOL      = 1.0e-6;     // accuracy for convergence test (minimum)
StepTOL  = 1.0e-6;     // accuracy for convergence test - size of the step 
XTOL     = 1.0e-6;     // accuracy for convergence test - improvement on x between two iterations
SubTOL   = 1.0e-6;     // accuracy for convergence test (line search)
Log      = %F;          // print info from line step, 1 = on / 0 = off
Spectral = %F;
h        = 4.0; //0.01;       // initial step length in line search
// Parameters specific to LBFGS
n        = 20;

TestCG    = %F;
TestBFGS  = %F;
TestLBFGS = %T;

Plot = %T;

deff('y=f(x)','y='+funcname+'(x,n);');
deff('y=df(x)','y=derivative(f,x,h=1e-2)'';');

x0 = zeros(n,1);

// lsm  = ls_secant;
lsm  = ls_dicho;
// lsm  = ls_goldsect;
// lsm  = ls_backtrack;
// lsm  = ls_polynom;

printf('\n*************** start example %s ***************\n\n', funcname);

lines(0);

/////////////////////////
// Limited memory BFGS //
/////////////////////////
if TestCG then
  printf('lbfgsdemo: CG - starting point: f(x) = %f - ||df(x)|| = %f\n', f(x0), norm(df(x0)));
  [x_opt, x_history] = optim_cg(x0, f, df, h, Log, ItMX, ls_ItMX, lsm, cgm, TOL, StepTOL, XTOL, SubTOL, Restart, Spectral);
  printf('lbfgsdemo: CG - solution found: f(x) = %f - ||df(x)|| = %f\n', f(x_opt), norm(df(x_opt)));

  NbIter = 0;
  for i=1:length(x_history)
    NbIter = NbIter + length(x_history(i));
  end
  printf('lbfgsdemo: CG - Number of call to the objective function = %d\n',NbIter);
end

if TestBFGS then
  printf('lbfgsdemo: BFGS - starting point: f(x) = %f - ||df(x)|| = %f\n', f(x0), norm(df(x0)));
  [x_opt, x_history] = optim_bfgs(x0, f, df, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL);
  printf('lbfgsdemo: BFGS - solution found: f(x) = %f - ||df(x)|| = %f\n', f(x_opt), norm(df(x_opt)));

  NbIter = 0;
  for i=1:length(x_history)
    NbIter = NbIter + length(x_history(i));
  end
  printf('lbfgsdemo: BFGS - Number of call to the objective function = %d\n',NbIter);
end

if TestLBFGS then
  if Plot then
    scf();
    xtitle('Comparaison','Iteration','F');
  end
  
  Legends(1) = 'm = 5';
  Legends(2) = 'm = 10';
  Legends(3) = 'm = 15';
  Legends(4) = 'm = 19';
  
  m     = [5,10,15,19];
  Style = ['r-','g-','k-','b-'];
  for i=1:length(m)
    printf('lbfgsdemo: LBFGS - starting point: f(x) = %f - ||df(x)|| = %f\n', f(x0), norm(df(x0)));
    [x_opt, x_history] = optim_lbfgs(x0, f, df, h, Log, ItMX, ls_ItMX, lsm, TOL, StepTOL, XTOL, SubTOL, Restart, m(i));
    printf('lbfgsdemo: LBFGS - solution found: f(x) = %f - ||df(x)|| = %f\n', f(x_opt), norm(df(x_opt)));

    NbIter = 0;
    for j=1:length(x_history)
      NbIter = NbIter + length(x_history(j));
    end
    printf('lbfgsdemo: LBFGS - Number of call to the objective function = %d\n',NbIter);
    if Plot then
      Index = 1;
      F = []; T = [];
      for j=1:length(x_history)
        F(Index) = f(x_history(j)(1));
        T(Index ) = Index;
        Index = Index + 1;
      end
      plot(T,F,Style(i));
    end
  end
  
  if Plot then
    legend(Legends);
  end
end



function [fmin, xmin, out] = cmaes( ... // output along Scilab function optim
    fitfun, x0, sigma, ...              // mandatory arguments
    stopfitness, stopfunevals, ...      // optional arguments
    mu, lambda, ...                     // population sizing
    displaymod, logmod, savemod, plotmod) // verbosity
// CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
// nonlinear function minimization. To be used under the terms of the
// GNU General Public License (http://www.gnu.org/copyleft/gpl.html).
//
// Author: Nikolaus Hansen, 2003. 
// e-mail: hansen[at]bionik.tu-berlin.de
// URL: http://www.bionik.tu-berlin.de/user/niko
// References: See end of file. Last change: October, 27, 2004

// Modification for scilab:
// Author: Yann COLLETTE, 2006. 
// e-mail: yann[dot]colletet[at]renault[dot]com
// URL: http://ycollette.free.fr

// Many more modifications and extensions in Scilab by Nikolaus Hansen, 2006

// Output :
// fmin : the objective function value of the best solution found
// xmin : best solution found 
// out: struct with all relevant data
// Input :
// f : the objective function (mandatory)
// x0 : the value of the initial point (mandatory)
// sigma : coordinate wise standard deviation (step size) (mandatory)
// stopfitness : stop if fitness < stopfitness (minimization) (optional)
// stopfunevals : stop after stopfunevals number of function evaluations (optional)
// lambda : size of the population (optional)
// mu : number of parents / point for recombination (optional)
// displaymod >= 1 : display some text every displaymod iteration
// logmod >= 0: log output data struct out every logmod iteration
// savemod >= 0: save everything every savemod iteration
// savemod >= 0: save everything every savemod iteration
// plotmod >= 0: plots everything every plotmod iteration
// with [plotmod] = resume(newvalue) after Ctrl-C the variable can be 
//   modified online. 

[nargout, nargin] = argn();

// --------------------  Initialization --------------------------------  

// mandatory arguments
if (~isdef('fitfun','local')) 
  error('CMAES: fitfun must be defined');
end
if (~isdef('x0','local')) 
  error('CMAES: x0 must be defined');
else
  N = size(x0,1); // number of objective variables/problem dimension
  xmean = x0; // objective variables initial point
end
if size(x0, 2) > 1
  error('CMAES: x0 must be a column vector');
end
if (~isdef('sigma','local')) 
  error('sigma must be defined'); // coordinate wise standard deviation (step size)
end

// optional arguments
if (~isdef('displaymod','local')) 
  displaymod = 100;
elseif displaymod
  printf('displaymod=%d\n', displaymod); 
end
if (~isdef('logmod','local')) 
  logmod = 1;
elseif displaymod
  printf('logmod=%d\n', logmod); 
end

if (~isdef('savemod','local')) 
  savemod = max(10, logmod);
elseif displaymod
  printf('savemod=%d\n', savemod); 
end

if ~isdef('plotmod','local')
  plotmod = max(100, logmod);
elseif displaymod
  printf('plotmod=%d\n', plotmod); 
end

if (~isdef('stopfitness','local') | isempty(stopfitness)) 
  stopfitness = -%inf; // stop if fitness < stopfitness (minimization)
elseif displaymod
  printf('stopfitness=%f\n', stopfitness); 
end
if (~isdef('stopfunevals','local')) 
  stopfunevals = 1e3*N^2;      // stop after stopfunevals number of function evaluations
elseif displaymod
  printf('stopfunevals=%.0g\n', stopfunevals); 
end

// --- Strategy parameter setting: Selection  
if (~isdef('lambda','local') | isempty(lambda) | lambda == 0) 
  lambda = 4+floor(3*log(N)); // population size, offspring number
elseif displaymod
  printf('lambda=%d\n', lambda); 
end
if (~isdef('mu','local') | isempty(mu) | mu == 0) 
  mu = floor(lambda/2); // number of parents/points for recombination
elseif displaymod
  printf('mu=%d\n', mu); 
end

weights = log(mu+1)-log(1:mu)';           // muXone array for weighted recombination
// lambda=12; mu=3; weights = ones(mu,1); // uncomment for (3_I,12)-ES
weights = weights/sum(weights);           // normalize recombination weights array
mueff   = sum(weights)^2/sum(weights.^2); // variance-effective size of mu

// Strategy parameter setting: Adaptation
cc     = 4/(N+4);                         // time constant for cumulation for covariance matrix
cs     = (mueff+2)/(N+mueff+3);           // time constant for cumulation for sigma control
mucov  = mueff;                           // size of mu used for calculating learning rate ccov
ccov   = (1/mucov) * 2/(N+1.4)^2 + (1-1/mucov) * ... // learning rate for covariance matrix
        ((2*mueff-1)/((N+2)^2+2*mueff));             //   
damps  = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; // damping for sigma 
                                                     // usually close to 1

// Initialize dynamic (internal) strategy parameters and constants
pc   = zeros(N,1); ps = zeros(N,1);  // evolution paths for C and sigma
B    = eye(N,N);                     // B defines the coordinate system
D    = eye(N,N);                     // diagonal matrix D defines the scaling
C    = B*D*(B*D)';                   // covariance matrix
chiN = N^0.5*(1-1/(4*N)+1/(21*N^2)); // expectation of ||N(0,I)|| == norm(randn(N,1))
                                     // exact: sqrt(2) * gamma((N+1)/2) / gamma(N/2) 
out.seed = rand('seed'); 

// -------------------- Generation Loop --------------------------------

if (displaymod) 
  printf('(%d/%d_W,%d)-CMAES \nEvals :  Function Value | Axis Ratio\n', mu, mu, lambda);
end

countiter = 0; 
counteval = 0;  // the next 40 lines contain the 20 lines of interesting code 
while counteval < stopfunevals
  countiter = countiter + 1;

  // Generate and evaluate lambda offspring
  arz = rand(N, lambda, 'normal');  // array of normally distributed mutation vectors
  for k=1:lambda,
    arx(:,k) = xmean + sigma * (B*D * arz(:,k));   // add mutation  // Eq. (1)
    if type(fitfun) == 13 | type(fitfun) == 11 // is a function
      arfitness(k) = fitfun(arx(:,k)); // objective function call
    elseif type(fitfun) == 15  | type(fitfun) == 16 // is a list
      error('not yet implemented, see function optim as guide');
      arfitness(k) = evalstr(fitfun(1), fitfun(2:$));
    else
      error('fitfun has invalid type');
    end
    counteval    = counteval+1;
  end

  // Sort by fitness and compute weighted mean into xmean
  [arfitness, arindex] = sort(-arfitness); // minimization
  arfitness = - arfitness;

  xmean = arx(:,arindex(1:mu))*weights;   // recombination, new mean value
  zmean = arz(:,arindex(1:mu))*weights;   // == sigma^-1*D^-1*B'*(xmean-xold)

  // Cumulation: Update evolution paths
  ps   = (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B * zmean);            // Eq. (4)
  hsig = norm(ps) / sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.5 + 1/(N+1);
  pc   = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (B * D * zmean); // Eq. (2)

  // Adapt covariance matrix C
  C = (1-ccov) * C ...                    // regard old matrix       // Eq. (3)
       + ccov * (1/mucov) * (pc*pc' ...   // plus rank one update
                             + (1-hsig) * cc*(2-cc) * C) ...
       + ccov * (1-1/mucov) ...           // plus rank mu update 
         * (B*D*arz(:,arindex(1:mu))) ...
         *  diag(weights) * (B*D*arz(:,arindex(1:mu)))';               

  // Adapt step size sigma
  sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));               // Eq. (5)
  
  // Update B and D from C
  // This is O(N^3). When strategy internal CPU-time is critical, the
  // next three lines can be executed only every (alpha/ccov/N)-th
  // iteration step, where alpha is e.g. between 0.1 and 10 
  C     = triu(C)+triu(C,1)';  // enforce symmetry
  [B,D] = spec(C);             // eigen decomposition, B==normalized eigenvectors
  D     = diag(sqrt(diag(D))); // D contains standard deviations now

  // end implementation of core algorithm

  // Set stopflag
  out.stopflag = list(); // for first iteration, otherwise it is empty anyway
  if max(diag(D)) / min(diag(D)) > 1e8
    out.stopflag($+1) = 'conditioncov';
  end
  if (arfitness(1) <= stopfitness) 
    out.stopflag($+1) = 'fitness'; 
  end
  if counteval >= stopfunevals
      out.stopflag($+1) = 'maxfunevals'; 
  end

  // Collect output data
  if logmod & (countiter < 10 ...
               | modulo(countiter-1,logmod) == 0 ...
               | ~isempty(out.stopflag)) 
      fmean = %nan; // clearly a hack, by default no additional fevals
      // Keep overall best solution
      if countiter == 1
	out.solutions.bestever.x = xmean;
	out.solutions.bestever.f = %inf;  // first hack
	out.solutions.bestever.evals = counteval;
	bestever = out.solutions.bestever;
      end
      out.evals = counteval;
      out.solutions.evals = counteval;
      out.solutions.mean.x = xmean;
      out.solutions.mean.f = fmean;
      out.solutions.mean.evals = counteval;
      out.solutions.recentbest.x = arx(:, arindex(1));
      out.solutions.recentbest.f = arfitness(1);
      out.solutions.recentbest.evals = counteval + arindex(1) - lambda;
      out.solutions.recentworst.x = arx(:, arindex($));
      out.solutions.recentworst.f = arfitness($);
      out.solutions.recentworst.evals = counteval + arindex($) - lambda;
      if arfitness(1) < out.solutions.bestever.f
	out.solutions.bestever.x = arx(:, arindex(1));
	out.solutions.bestever.f = arfitness(1);
	out.solutions.bestever.evals = counteval + arindex(1) - lambda;
	bestever = out.solutions.bestever;
      end
      
      out.evo.evals($+1) = counteval;
      out.evo.iterations($+1) = countiter;
      out.evo.mean.f($+1) = fmean;
      out.evo.mean.evals($+1) = counteval;
      out.evo.recentbest.f($+1) = arfitness(1); // reevaluations would make an array
      out.evo.recentbest.evals($+1) = counteval + arindex(1) - lambda;
      out.evo.recentworst.f($+1) = arfitness($); 
      out.evo.recentworst.evals($+1) = counteval + arindex($) - lambda;
      if countiter == 1
	out.evo.mean.x(:,1) = xmean;
	out.evo.recentbest.x(:,1) = arx(:,arindex(1));
	out.evo.recentworst.x(:,1) = arx(:,arindex($));
      else
	out.evo.mean.x(:,$+1) = xmean;
	out.evo.recentbest.x(:,$+1) = arx(:,arindex(1));
	out.evo.recentworst.x(:,$+1) = arx(:,arindex($));
      end	
    
      // Single Parameters
      out.evo.param.evals($+1) = counteval;
      out.evo.param.iterations($+1) = countiter;
      out.evo.param.sigma($+1) = sigma;
      out.evo.param.maxstd($+1) = sigma * sqrt(max(diag(C)));
      out.evo.param.minstd($+1) = sigma * sqrt(min(diag(C)));
      [ignore out.evo.param.maxstdidx($+1)] = max(diag(C));
      [ignore out.evo.param.minstdidx($+1)] = min(diag(C));
      out.evo.param.maxD($+1) = max(diag(D));
      out.evo.param.minD($+1) = min(diag(D));
      out.evo.param.comment = '';
      
      // Parameter Arrays
      out.evoParamArr.evals($+1) = counteval;
      out.evoParamArr.iterations($+1) = countiter;
      out.evoParamArr.sigma($+1) = sigma; // for convenience and completeness
      if countiter == 1
	out.evoParamArr.diagD(:,1) = diag(D); 
	out.evoParamArr.stds(:,1) = sigma * sqrt(diag(C));
	out.evoParamArr.Bmax(:,1) = B(:,out.evo.param.maxstdidx($)); 
	out.evoParamArr.Bmin(:,1) = B(:,out.evo.param.minstdidx($)); 
      else
	out.evoParamArr.diagD(:,$+1) = diag(D); 
	out.evoParamArr.stds(:,$+1) = sigma * sqrt(diag(C));
	out.evoParamArr.Bmax(:,$+1) = B(:,out.evo.param.maxstdidx($)); 
	out.evoParamArr.Bmin(:,$+1) = B(:,out.evo.param.minstdidx($)); 
      end
      out.evoParamArr.comment = '';
  end	

  // save output data
  if (savemod & (modulo(countiter-1,savemod) == 0) | ~isempty(out.stopflag)) 
    save('variablescmaes.dat');
  end

  // plot output data
  if plotmod & (countiter < 4 ...
                | modulo(countiter-1,plotmod) == 0 ...
                | ~isempty(out.stopflag)) 
    plotcmaes(out); 
  end

  // print / display 
  if (displaymod & (~isempty(out.stopflag) | modulo(countiter-1,displaymod) == 0)) 
    printf('%5d : %15.8e | %.2e\n', counteval, arfitness(1), max(diag(D))/min(diag(D)));
  end

  if ~isempty(out.stopflag)
    break;
  end

end // while, end generation loop

// -------------------- Final Assignments ---------------------------------

xmin = arx(:, arindex(1)); // Return best point of last generation.
                           // Notice that xmean is expected to be even
                           // better.
fmin = arfitness(arindex(1));

if displaymod
  for i = 1:length(out.stopflag)
    disp(strcat(['stop on ', out.stopflag(i)]));
  end
end
endfunction

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plotcmaes(out)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// input: out-struct from function cmaes

drawlater; 
clf; 

//////////////////
subplot(2,2,1); 
xgrid;
fmin = min(out.evo.recentbest.f);
dat.logy = log10(out.evo.recentbest.f - fmin + 1e-49);
idx = out.evo.recentbest.f == fmin;
fmin2 = min(dat.logy(~idx));
xtitle("Function Values minus f_min, Sigma", "", "log10(value)");
xtitle("Function Values minus f_min, Sigma", ...
       'f_recent='+string(out.evo.recentbest.f($)));

// plot marker(s) for best function value
plot(out.evo.recentbest.evals(idx), log10(fmin + 1e-49), 'r*');
legend('f_min=' + string(fmin));
plot(out.evo.recentbest.evals(idx), fmin2 - 0.1*(max(dat.logy)-fmin2), 'b*');

// plot abs function value in red
plot2d(out.evo.recentworst.evals, log10(abs(out.evo.recentworst.f) + 1e-49));

// plot sigma (green)
plot(out.evo.param.evals, log10(out.evo.param.sigma), 'g');
plot(out.evo.recentbest.evals, log10(1e-49 + abs(out.evo.recentbest.f)), 'r');
  // legend([strcat(['f_min=', string(fmin)]), 'sigma']);

// plot function value differences disregarding all fmins
istart = 1;
for i = find(idx)
  plot2d(out.evo.recentbest.evals(istart:i-1), ...
         dat.logy(istart:i-1), style=2); // blue style
  istart = i+1;
end
plot2d(out.evo.recentbest.evals(istart:$), ...
       dat.logy(istart:$), style=2); // blue style

//////////////////
subplot(2,2,2); 
xgrid;
xtitle(sprintf("Object Variables (%d-D)", size(out.solutions.bestever.x,1)));
plot(out.evo.mean.evals, out.evo.mean.x);

//////////////////
subplot(2,2,3);
xtitle("Standard Deviations", 'function evaluations', 'log10(value)');
plot(out.evoParamArr.evals, log10(out.evoParamArr.stds+1e-49));
xgrid;

//////////////////
subplot(2,2,4);
xtitle("Principle Axis Lengths", 'function evaluations', 'log10(value)');
plot(out.evoParamArr.evals, log10(out.evoParamArr.diagD+1e-49));
xgrid;
drawnow; 

endfunction 

// ---------------------------------------------------------------  
////// REFERENCES
//
// The equation numbers refer to 
// Hansen, N. and S. Kern (2004). Evaluating the CMA Evolution
// Strategy on Multimodal Test Functions.  Eighth International
// Conference on Parallel Problem Solving from Nature PPSN VIII,
// Proceedings, pp. 282-291, Berlin: Springer. 
// (http://www.bionik.tu-berlin.de/user/niko/ppsn2004hansenkern.pdf)
// 
// Further references:
// Hansen, N. and A. Ostermeier (2001). Completely Derandomized
// Self-Adaptation in Evolution Strategies. Evolutionary Computation,
// 9(2), pp. 159-195.
// (http://www.bionik.tu-berlin.de/user/niko/cmaartic.pdf).
//
// Hansen, N., S.D. Mueller and P. Koumoutsakos (2003). Reducing the
// Time Complexity of the Derandomized Evolution Strategy with
// Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
// 11(1).  (http://mitpress.mit.edu/journals/pdf/evco_11_1_1_0.pdf).
//


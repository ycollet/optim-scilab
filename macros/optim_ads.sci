function [x, x_history] = optim_ads(f, x, Tol, ItMX, Log, P)
//optim_ads  Alternating directions method for direct search optimization.
//        [x, x_history] = optim_ads(f, x, Tol, ItMX, P) attempts to
//        maximize the function f, using the starting vector x0.
//        The alternating directions direct search method is used.
//        Output arguments:
//               x    = vector yielding largest function value found,
//               x_history = list of each computed simplex
//        The iteration is terminated when either
//               - the relative increase in function value between successive
//                 iterations is <= Tol (default 1e-3),
//               - ItMX function evaluations have been performed
//                 (default inf, i.e., no limit), or
//        Progress of the iteration is not shown if Log = %F (default %F).
//        By default, the search directions are the co-ordinate directions.
//        The columns of a fifth parameter matrix P specify alternative search
//        directions (P = EYE is the default).
//        NB: x0 can be a matrix. 

//     Reference:
//     N. J. Higham, Optimization by direct search in matrix computations,
//        SIAM J. Matrix Anal. Appl, 14(2): 317-333, 1993.
//     N. J. Higham, Accuracy and Stability of Numerical Algorithms,
//        Second edition, Society for Industrial and Applied Mathematics,
//        Philadelphia, PA, 2002; sec. 20.5.

[nargout,nargin] = argn();

x0 = x(:);  // Work with column vector internally.
n = length(x0);

mu = 1e-4;  // Initial percentage change in components.
nstep = 25; // Max number of times to double or decrease h.

x_history_defined = (nargout==2);

if x_history_defined then
  x_history = list();
end

// Set up convergence parameters.
if ~isdef('Tol','local') then
  Tol = 1e-3;
end
if ~isdef('ItMX','local') then
  ItMX = 1000;
end
if ~isdef('Log','local') then
  Log  = %F;
end
if ~isdef('P','local') then
   P = eye(n,n);             // Matrix of search directions.
else
   if (and(size(P)~=[n n])) then  // Check for common error.
      error('P must be of dimension the number of elements in x0.')
   end
end

fmin = f(x); nf = 1;
if Log then printf('f(x0) = %9.4e\n', fmin), end

steps = zeros(n,1);
it = 0; y = x0;

while 1    // Outer loop.
  it = it+1;
  if Log then printf('Iter %2.0f  (nf = %2.0f)\n', it, nf), end
  fmin_old = fmin;

  for i=1:n  // Loop over search directions.
      pi = P(:,i);
      flast = fmin;
      yi = y;
      h = sign(pi'*yi)*norm(pi.*yi)*mu;   // Initial step size.
      if h == 0 then h = max(norm(yi,%inf),1)*mu; end
      y = yi + h*pi;
      x(:) = y; fnew = f(x); nf = nf + 1;
      if fnew < fmin then
         fmin = fnew;
         h = 2*h; lim = nstep; k = 1;
      else
         h = -h; lim = nstep+1; k = 0;
      end

      for j=1:lim
        y = yi + h*pi;
        x(:) = y; fnew = f(x); nf = nf + 1;
        
        if x_history_defined then
          x_history($+1) = x;
        end
      
        if fnew >= fmin then break, end
        fmin = fnew; k = k + 1;
        h = 2*h;
      end

     steps(i) = k;
     y = yi + 0.5*h*pi;
     if k == 0 then y = yi; end

     if x_history_defined then
       x_history($+1) = x;
     end

     if Log then
        printf('Comp. = %2.0f,  steps = %2.0f,  f = %9.4e', i, k, fmin)
        printf('  (%2.1f%%)\n', 100*abs(fmin-flast)/(abs(flast)+%eps))
     end

     if nf >= ItMX then
        if Log then
           printf('Max no. of function evaluations exceeded...quitting\n')
        end
        x(:) = y;
        return;
     end
  end  // Loop over search directions.

  if isequal(steps,zeros(n,1))
     if Log then printf('Stagnated...quitting\n'), end
     x(:) = y; 
     return;
  end

  if (fmin_old-fmin) <= Tol*abs(fmin_old)
     if Log then printf('Function values ''converged''...quitting\n'), end
     x(:) = y;
     return;
  end

end //////////// Of outer loop.
endfunction

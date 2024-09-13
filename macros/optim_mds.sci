function [x, x_history] = optim_mds(fun, x, SimplexSize, ItMX, SimplexShape, Log)
//optim_mds  Multidirectional search method for direct search optimization (taken from mdsmax of octave).
//        [x, x_history] = optim_mds(fun, x, SimplexSize, ItMX, SimplexShape, Log) attempts to
//        minimize the function fun, using the starting vector x0.
//        The method of multidirectional search is used.
//        Output arguments:
//               x    = vector yielding largest function value found,
//               x_history = list of x values computed
//        The form of the initial simplex is determined by SimplexShape:
//          SimplexShape = %F: regular simplex (sides of equal length, the default),
//          SimplexShape = %T: right-angled simplex.
//        The size of the initial simplex is determined by SimplexSize (default value: 1e-3)
//        The maximum number of iterations is determined by ItMX (default value: 1000)
//        You can trace the behavior of optim_mds using Log (default value: %F - no trace)
//        NB: x0 can be a matrix.  

// This implementation uses 2n^2 elements of storage (two simplices), where x0
// is an n-vector.  It is based on the algorithm statement in [2, sec.3],
// modified so as to halve the storage (with a slight loss in readability).

// From Matrix Toolbox 
// Copyright (C) 2002 N.J.Higham
// www.maths.man.ac.uk/~higham/mctoolbox
// distributed under the terms of the GNU General Public License
//
// Modifications for octave by A.Adler 2003
// $Id: mdsmax.m,v 1.4 2005/05/25 03:43:41 pkienzle Exp $

// References:
// [1] V. J. Torczon, Multi-directional search: A direct search algorithm for
//     parallel machines, Ph.D. Thesis, Rice University, Houston, Texas, 1989.
// [2] V. J. Torczon, On the convergence of the multidirectional search
//     algorithm, SIAM J. Optimization, 1 (1991), pp. 123-145.
// [3] N. J. Higham, Optimization by direct search in matrix computations,
//     SIAM J. Matrix Anal. Appl, 14(2): 317-333, 1993.
// [4] N. J. Higham, Accuracy and Stability of Numerical Algorithms,
//        Second edition, Society for Industrial and Applied Mathematics,
//        Philadelphia, PA, 2002; sec. 20.5.

[nargout,nargin] = argn();

x_history_defined = (nargout==2);

x0 = x(:);  // Work with column vector internally.
n  = length(x0);

mu    = 2;   // Expansion factor.
theta = 0.5; // Contraction factor.

// Set up convergence parameters etc.
if ~isdef('SimplexSize','local') then
  SimplexSize = 1e-3;
end
if ~isdef('ItMX','local') then
  ItMX = 1000;
end
if ~isdef('SimplexShape','local') then
  SimplexShape = %F;
end
if ~isdef('Log','local') then
  Log = %F;
end

if x_history_defined then
  x_history = list();
end

V = [zeros(n,1) eye(n,n)]; T = V;
f = zeros(n+1,1); ft = f;
V(:,1) = x0; f(1) = -fun(x);
fmax_old = f(1);

if Log then printf('f(x0) = %9.4e\n', -f(1)); end

k = 0; m = 0;

// Set up initial simplex.
scale = max(norm(x0,%inf),1);
if SimplexShape == %F then
   // Regular simplex - all edges have same length.
   // Generated from construction given in reference [18, pp. 80-81] of [1].
   alpha = scale / (n*sqrt(2)) * [ sqrt(n+1)-1+n  sqrt(n+1)-1 ];
   V(:,2:n+1) = (x0 + alpha(2)*ones(n,1)) * ones(1,n);
   for j=2:n+1
       V(j-1,j) = x0(j-1) + alpha(1);
       x(:) = V(:,j); f(j) = -fun(x);
   end
else
   // Right-angled simplex based on co-ordinate axes.
   alpha = scale*ones(n+1,1);
   for j=2:n+1
       V(:,j) = x0 + alpha(j)*V(:,j);
       x(:) = V(:,j); f(j) = -fun(x);
   end
end

if x_history_defined then
  x_history($+1) = list();
  for i=1:size(V,2)
    x_history($)($+1) = V(:,i);
  end
end

nf = n+1;
_size = 0;         // Integer that keeps track of expansions/contractions.
flag_break = 0;   // Flag which becomes true when ready to quit outer loop.

while 1    //////////// Outer loop.
  k = k+1;

  // Find a new best vertex  x  and function value  fmax = f(x).
  [fmax,j] = max(f);
  V(:,[1 j]) = V(:,[j 1]); v1 = V(:,1);
  f([1 j]) = f([j 1]);
  if Log then
     printf('Iter. %2.0f,  inner = %2.0f,  size = %2.0f,  ', k, m, _size);
     printf('nf = %3.0f,  f = %9.4e  (%2.1f%%)\n', nf, -fmax, 100*(fmax-fmax_old)/(abs(fmax_old)+%eps));
  end
  fmax_old = fmax;

  m = 0;
  while 1   ////// Inner repeat loop.
      m = m+1;

      // Stopping Test 2 - too many f-evals?
      if nf >= ItMX then
         msg = ['Max no. of function evaluations exceeded...quitting\n'];
         flag_break = 1; 
         break  // Quit.
      end

      // Stopping Test 3 - converged?   This is test (4.3) in [1].
      size_simplex = norm(V(:,2:n+1)- v1(:,ones(1,n)),1) / max(1, norm(v1,1));
      if size_simplex <= SimplexSize then
         msg = sprintf('Simplex size %9.4e <= %9.4e...quitting\n', size_simplex, SimplexSize);
         flag_break = 1; 
         break  // Quit.
      end

      for j=2:n+1      // ---Rotation (reflection) step.
          T(:,j) = 2*v1 - V(:,j);
          x(:) = T(:,j); ft(j) = -fun(x);
      end
    
      nf = nf + n;

      replaced = ( max(ft(2:n+1)) > fmax );

      if replaced then
         for j=2:n+1   // ---Expansion step.
             V(:,j) = (1-mu)*v1 + mu*T(:,j);
             x(:) = V(:,j); f(j) = -fun(x);
         end
         nf = nf + n;
         // Accept expansion or rotation?
         if max(ft(2:n+1)) > max(f(2:n+1)) then
            V(:,2:n+1) = T(:,2:n+1);  f(2:n+1) = ft(2:n+1);  // Accept rotation.
         else
            _size = _size + 1;  // Accept expansion (f and V already set).
         end
      else
         for j=2:n+1   // ---Contraction step.
             V(:,j) = (1+theta)*v1 - theta*T(:,j);
             x(:) = V(:,j); f(j) = -fun(x);
         end
         nf = nf + n;
         replaced = ( max(f(2:n+1)) > fmax );
         // Accept contraction (f and V already set).
         _size = _size - 1;
      end

      if x_history_defined then
        x_history($+1) = list();
        for i=1:size(V,2)
          x_history($)($+1) = V(:,i);
        end
      end

      if replaced then break; end
      if Log & modulo(m,10) == 0 then printf('        ...inner = %2.0f...\n',m); end
      end ////// Of inner repeat loop.

  if flag_break then break; end
end //////////// Of outer loop.

// Finished.
if Log then printf(msg); end
x(:) = v1;
endfunction

function [x_opt, x_history] = optim_rosen(x0, rosen_f, h, Log, ItMX, TOL, StepTOL, XTOL, cgm, _alpha, _beta)
//
//////////////////
// Informations //
//////////////////
//
// Steepest descent optimization method
// 
// x0       : initial starting point (must be a column vector)
// rosen_f  : objective function
// h        : size of the initial step (optional parameter)
// Log      : verbose mode or not (%T / %F) (optional parameter)
// ItMX     : number of steepest descent step (optional parameter)
// ls_ItMX  : number of iterations of the line search (optional parameter)
// TOL      : accuracy for convergence test - derivatives (optional parameter)
// XTOL     : accuracy for convergence test - improvement on x between two iterations (optional parameter)
// cgm      : type of restart basis generation: rosen0, rosen1, rosen2
// _alpha   : Increase step size coefficient
// _beta    : Decrease step size coefficient
//

///////////////////////////
// Testing the arguments //
///////////////////////////

[nargout, nargin] = argn();

x_history_defined = (nargout==2);

if ~isdef('x0','local') then
  error('optim_rosen: x0 is mandatory');
end

if (~isdef('h','local')) then
  h = 0.00125; // initial step length in line search
end

if (~isdef('Log','local')) then
  Log = %F; // print info from line step, %T = on / %F = off
end

if (~isdef('ItMX','local')) then
  ItMX = 50; // maximum number of descent steps
end

if (~isdef('TOL','local')) then
  TOL = 1.0e-4; // accuracy for convergence test - derivatives 
end

if (~isdef('XTOL','local')) then
  XTOL = 1.0e-4; // accuracy for convergence test - improvement on x between two iterations
end

if (~isdef('cgm','local')) then
  cgm = 'rosen1';
end

if (~isdef('_alpha','local')) then
  _alpha = 2;
end

if (~isdef('_beta','local')) then
  _beta = 0.5;
end

if ~isdef('rosen_f','local') then
  error('optim_rosen: rosen_f is mandatory');
else
  if typeof(rosen_f)=='list' then
    deff('y=_rosen_f(x)','y=rosen_f(1)(x, rosen_f(2:$))');
  else
    deff('y=_rosen_f(x)','y=rosen_f(x)');
  end
end

////////////////// Rosenbrock method

n      = size(x0,1);
bigbnd = 1e10;
It     = 0;
dirchg = 0;
xi     = eye(n,n);
A      = eye(n,n);
x      = zeros(n,1);
x_test = zeros(n,1);
x1     = zeros(n,1);
epsilon = 1e-5;

if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end

// Copy of the initial point
x1 = x0;

// Initialization of the basis size steps
d      = h*ones(n,1);
lambda = zeros(n,1);

// Computation of the first objective function value

yfirstfirst = _rosen_f(x1);
yfirst      = yfirstfirst;
It          = It + 1; 

if (x_history_defined) then
  x_history($+1) = x1;
end

// We set restart at true to pass the first test
FirstTime = %T;
Restart   = %T;

while (((Restart) & (It<ItMX)) | (FirstTime)) do
  ybest = yfirstfirst;
  
  while ((ybest<yfirst) | (FirstTime)) do
    FirstTime = %F;
    yfirst = ybest;

    // We test each sens
    for l=1:2
      // we test each directions
      for i=1:n
        // Construction of the vector displacement using the xi basis and the d coordinates
        for j=1:n
          x_test(j) = x1(j) + d(j)*xi(j,i);
        end
            
        // Test of the direction
        ycurrent = _rosen_f(x_test);
        It       = It + 1;
        
        if (x_history_defined) then
          x_history($+1) = x_test;
        end

        if (ycurrent < ybest) then
          // If the vector displacement reduce the value of the objective function,
          // we store xcurrent as the new x1 and increase the size of d
          // lambda holds the total displacement (displacement from x0 to xopt)
          lambda(i) = lambda(i) + d(i);
          d(i)      = d(i) * _alpha;
          ybest     = ycurrent;
          x1        = x_test;
        else
          // If it doesn't reduce the value of the objective function, we search
          // in the other direction and reduce the size of d
          d(i) = d(i) * (- _beta);
        end // if else
        
        if (It>ItMX) then
          break;
        end
      end // for
      
      if (It>ItMX) then
        break;
      end
    end // for
  end // while
    
  if (It>ItMX) then
    break;
  end

  _mini = bigbnd;
  // We check that the size step is not too small
  _mini = min([_mini abs(d')]);
  // If the size is too small, we restart the search using canonical basis
  Restart = (_mini > epsilon);
    
  if (ybest < yfirstfirst) then
    _mini = bigbnd;

    // If the displacement between x1 and x is too small, we order a restart
    _mini = min([_mini abs(x1' - x')]);
    Restart = (Restart | (_mini > epsilon));
    
    if (Restart) then
      printf('optim_rosen: restart\n');
      // nous avons:
      // x1[j]-x[j]=(somme sur i de) lambda[i]*xi[i][j];
      dirchg = dirchg + 1; // count of the number of direction change
      
	    if (cgm=='rosen0') then  
        // Second, using the search results, we build another vector base
        // Building the A set of directions
        for i=1:n
          for j=i:n
            A(:,i) = A(:,i) + d(i) * xi(:,i);
          end
        end
        // Building the new basis
        B = A;
        for i=1:n
          Aux = A(:,i);
          if (i>1) then
            for j=1:i-1
              Aux = Aux - A(:,j)'*xi(:,j)*xi(:,j);
            end
          end
          B(:,i) = Aux / norm(Aux);    
        end
        A = B;
      end
	    
      if (cgm=='rosen1') then  
        // Second, using the search results, we build another vector base
        // Building the A set of directions using Swann modification
        // We sort all the steps by decreasing order
        [s, k] = sort(d);
        A = A(:,k);
        d = d(k);
      
        for i=1:n
          for j=i:n
            if (d(i)>StepTOL) then
              A(:,i) = A(:,i) + d(i) * xi(:,i);
            else
              // If the step is too small, we keep the old direction
              A(:,i) = xi(:,i);
            end
          end
        end
        // Building the new basis
        B = A;
        for i=1:n
          Aux = A(:,i);
          if (i>1) then
            for j=1:i-1
              Aux = Aux - A(:,j)'*xi(:,j)*xi(:,j);
            end
          end
          B(:,i) = Aux / norm(Aux);    
        end
        A = B;
      end
	    
	    if (cgm=='rosen2') then  
        A(:,n) = lambda(n)*xi(:,n);
      
        for k=n-1:-1:1
          A(:,k) = A(:,k+1) + lambda(k)*xi(:,k);
        end // for
	      
  	    t(n) = lambda(n)*lambda(n);
	      for i=n-1:-1:1
          t(i) = t(i+1) + lambda(i)*lambda(i);
        end // for
  	    for i=n:-1:2
          div = sqrt(t(i-1)*t(i));
          if (div~=0) then
            xi(:,i) = (lambda(i-1)*A(:,i) - xi(:,i-1)*t(i))/div;
          end // If
        end // for
        div = sqrt(t(1));
        if (div~=0) then
          xi(:,1) = A(:,1)/div;
        end
      end // If
    end  // if
  else
    // There was no improvement during the preceding search
    y = _rosen_f(x0);
    
    if (x_history_defined) then
      x_history($+1) = x0;
    end

    printf('\noptim_rosen: no improvement during iteration %d - after %d direction change:\n', It, dirchg);
    printf('%10s\n', 'f(x)');
    printf('%10.2e\n', y);
    FirstTime = %T;
  end // if

  Restart = %T;
  
  x0 = x1;
    
  if (norm(fs(x0))<TOL) then
    y = _rosen_f(x0);

    if (x_history_defined) then
      x_history($+1) = x0;
    end

    printf('\noptim_rosen: minimum found in iteration %d - after %d direction change :\n', It, dirchg);
    printf('%10s\n', 'f(x)');
    printf('%10.2e\n', y);
    x_opt = x0;
    return;
  else
    y = _rosen_f(x0);

    if (x_history_defined) then
      x_history($+1) = x0;
    end

    printf('optim_rosen: %10s\n', 'f(x)');
    printf('%10.2e\n', y);
  end // if else
end // while

x_opt = x0;
endfunction


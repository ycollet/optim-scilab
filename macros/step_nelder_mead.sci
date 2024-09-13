function [x_next, data_next, eval_Func, f_hist, x_hist] = step_nelder_mead(f_current, x_current, data_current, nm_mode, Log, kelley_restart, kelley_alpha)

[nargout, nargin] = argn();

ItMX = 100;
Tol  = 0.0;
n    = size(x_current,1);

if (~isdef('Log','local')) then
  Log = %F;
end
if (~isdef('kelley_restart','local')) then
  kelley_restart = %F;
end
if (~isdef('kelley_alpha','local')) then
  kelley_alpha = 1e-4;
end

f_hist_is_defined = (nargout>=4);
x_hist_is_defined = (nargout==5);

Compute = %F;

while ~Compute
  if (nm_mode=='init') then
    n = size(x_current,1);
    if (size(x_current,2)~=n+1) then
      printf('nelder_mead: you must give %d starting points\n',n+1);
      abort;
    end
    // Initialisation of the state structure
    data_next = list();
    data_next(25) = list();
    data_next(26) = list();
    for i=1:n+1
      data_next(25)($+1) = f_current(i);
      data_next(26)($+1) = x_current(:,i);
    end
    data_next(1)  = f_current(1);
    data_next(2)  = f_current(2);
    data_next(3)  = f_current(3);
    data_next(4)  = x_current(:,1);
    data_next(5)  = x_current(:,2);
    data_next(6)  = x_current(:,3);
    data_next(7)  = 0;  // Case
    data_next(8)  = -1;  // f_Reflexion
    data_next(9)  = []; // x_Reflexion
    data_next(10) = -1;  // f_Expansion
    data_next(11) = []; // x_Expansion
    data_next(12) = -1;  // f_C1
    data_next(13) = []; // x_C1
    data_next(14) = -1;  // f_C2;
    data_next(15) = []; // x_C2
    data_next(16) = -1;  // f_Contr;
    data_next(17) = []; // x_Contr;
    data_next(18) = -1;  // f_Shrink
    data_next(19) = []; // x_Shrink
    data_next(20) = -1;  // f_Median
    data_next(21) = []; // x_Median
    data_next(22) = list(); // f_hist
    data_next(23) = list(); // x_hist
    data_next(24) = 0; // eval_Func

    eval_Func = 0;
    f_nm = [];
    for i=1:n+1
      f_nm(i) = f_current(i);
    end

    [f_sort, I] = sort(f_nm);

    f_Worth = f_nm(I(1));
    f_Good  = f_nm(I($-1));
    f_Best  = f_nm(I($));
    x_Worth = data_next(26)(I(1));
    x_Good  = data_next(26)(I($-1));
    x_Best  = data_next(26)(I($));

    data_next(1) = f_Worth;
    data_next(2) = f_Good;
    data_next(3) = f_Best;
    data_next(4) = x_Worth;
    data_next(5) = x_Good;
    data_next(6) = x_Best;

    if (x_hist_is_defined) then
      data_next(23)($+1)  = list();
      for i=1:n+1
        data_next(23)($)($+1) = data_next(26)(i);
      end
    end

    if (f_hist_is_defined) then
      data_next(22)($+1)  = list();
      for i=1:n+1
        data_next(22)($)($+1) = data_next(25)(i);
      end
    end

    data_current = list();
    data_current(7) = 0; // Fake variable to allow the first computation after initialization

    f_Reflexion = -1;
    x_Reflexion = [];
    f_Expansion = -1;
    x_Expansion = [];
    f_C1        = -1;
    x_C1        = [];
    f_C2        = -1;
    x_C2        = [];
    f_Contr     = -1;
    x_Contr     = [];
    f_Shrink    = -1;
    x_Shrink    = [];
    f_Median    = -1;
    x_Median    = [];

    Case        = 1; // Starting position of the algorithm
  end

  if (nm_mode=='run') then
    data_next = list();

    f_Worth = data_current(1);
    f_Good  = data_current(2);
    f_Best  = data_current(3);
    x_Worth = data_current(4);
    x_Good  = data_current(5);
    x_Best  = data_current(6);

    f_Reflexion = data_current(8);
    x_Reflexion = data_current(9);
    f_Expansion = data_current(10);
    x_Expansion = data_current(11);
    f_C1        = data_current(12);
    x_C1        = data_current(13);
    f_C2        = data_current(14);
    x_C2        = data_current(15);
    f_Contr     = data_current(16);
    x_Contr     = data_current(17);
    f_Shrink    = data_current(18);
    x_Shrink    = data_current(19);
    f_Median    = data_current(20);
    x_Median    = data_current(21);

    eval_Func   = data_current(24);

    Case        = data_current(7);

    data_next   = data_current;
  end

  if (nm_mode=='exit') then
    x_next    = data_current(6);
    data_next = data_current(3);
    f_hist    = data_current(22);
    x_hist    = data_current(23);
    return;
  end

  if (Log) then
    printf('step_nelder_mead: initialization\n');
    printf('f_Worth = %f f_Good = %f f_Best = %f\n', f_Worth, f_Good, f_Best);
  end

  select(Case)
  case 1 then
    // Calcul de R
    x_Median      = (x_Best + x_Good) / 2.0
    x_Reflexion   = 2.0 * x_Median - x_Worth;
    x_next        = x_Reflexion;
    data_next(9)  = x_Reflexion;
    data_next(21) = x_Median;
    data_next(7)  = 2;
    Compute       = %T;
  case 2 then
    // On récupère f(R)
    f_Reflexion   = f_current;
    data_next(8)  = f_Reflexion;
    data_next(7)  = 3; // Next Case
    eval_Func     = eval_Func + 1;
    data_next(24) = eval_Func;
    Compute       = %F;
  case 3 then
    if (f_Reflexion<f_Good) then
      data_next(7) = 4; // Next Case
    else
      data_next(7) = 11; // Next Case
    end
    Compute = %F;
  case 4 then
    if (f_Best<f_Reflexion) then
      data_next(7) = 5; // Next Case
    else
      data_next(7) = 6;
    end
    Compute = %F;
  case 5 then
    f_Worth       = f_Reflexion;
    x_Worth       = x_Reflexion;
    data_next(1)  = f_Worth;
    data_next(4)  = x_Worth;
    data_next(7)  = 22; // Next Case (we start again), but before, we rank the points
    Compute = %F;
  case 6 then
    // Computation of E
    x_Expansion   = 2.0 * x_Reflexion - x_Median;
    x_next        = x_Expansion;
    data_next(11) = x_Expansion;
    data_next(7)  = 7;
    Compute       = %T;
  case 7 then
    // Computation of f(E)
    f_Expansion   = f_current;
    data_next(10) = f_Expansion;
    data_next(7)  = 8;
    eval_Func     = eval_Func + 1;
    data_next(24) = eval_Func;
    Compute       = %F;
  case 8 then
    if (f_Expansion<f_Best) then
      data_next(7) = 9;
    else
      data_next(7) = 10;
    end
    Compute = %F;
  case 9 then
    // We replace W by E
    x_Worth      = x_Expansion;
    f_Worth      = f_Expansion;
    data_next(1) = f_Worth;
    data_next(4) = x_Worth;
    data_next(7) = 22; // Next Case: we start again but before, we rank the points
    Compute      = %F;
  case 10 then
    // We replace W by R
    x_Worth       = x_Reflexion;
    f_Worth       = f_Reflexion;
    data_next(1)  = f_Worth;
    data_next(4)  = x_Worth;
    data_next(7)  = 22; // Next Case: we start again but before, we rank the points
    Compute       = %F;
  case 11 then
    // if f(R) < f(W)
    if (f_Reflexion<f_Worth) then
      data_next(7) = 12;
    else
      data_next(7) = 13;
    end
    Compute = %F;
  case 12 then
    // We replace W by R
    x_Worth      = x_Reflexion;
    f_Worth      = f_Reflexion;
    data_next(1) = f_Worth;
    data_next(4) = x_Worth;
    data_next(7) = 13; // Next Case: we start again
    Compute      = %F;
  case 13 then
    // Computation of C1
    x_C1          = (x_Median + x_Worth) / 2.0;
    x_next        = x_C1;
    data_next(13) = x_C1;
    data_next(7)  = 14;
    Compute = %T;
  case 14 then
    // Computation of f(C1)
    f_C1          = f_current;
    data_next(12) = f_C1;
    data_next(7)  = 15;
    eval_Func     = eval_Func + 1;
    data_next(24) = eval_Func;
    Compute       = %F;
  case 15 then
    // Computation of C2
    x_C2          = (x_Median + x_Reflexion) / 2.0;
    x_next        = x_C2;
    data_next(15) = x_C2;
    data_next(7)  = 16;
    Compute       = %T;
  case 16 then
    // Computation of f(C2)
    f_C2          = f_current;
    data_next(14) = f_C2;
    data_next(7)  = 17;
    eval_Func     = eval_Func + 1;
    data_next(24) = eval_Func;
    Compute       = %F;
  case 17 then
    if (f_C1<f_C2) then
      f_Contr = f_C1;
      x_Contr = x_C1;
    else
      f_Contr = f_C2;
      x_Contr = x_C2;
    end
    data_next(16) = f_Contr;
    data_next(17) = x_Contr;
    if (f_Contr<f_Worth) then
      data_next(7) = 21;
    else
      data_next(7) = 18;
    end
    Compute = %F;
  case 18 then
    // Computation of S
    x_Shrink      = (x_Best + x_Worth) / 2.0;
    x_next        = x_Shrink;
    data_next(19) = x_Shrink;
    data_next(7)  = 19;
    Compute       = %T;
  case 19 then
    // Computation of f(S)
    f_Shrink      = f_current;
    data_next(18) = f_Shrink;
    data_next(7)  = 20;
    eval_Func     = eval_Func + 1;
    data_next(24) = eval_Func;
    Compute       = %F;
  case 20 then
    // We replace W by S and G by M
    x_Worth      = x_Shrink;
    f_Worth      = f_Shrink;
    x_Good       = x_Median;
    f_Good       = f_Median;
    data_next(2) = f_Good;
    data_next(1) = f_Worth;
    data_next(5) = x_Good;
    data_next(4) = x_Worth;
    data_next(7) = 22; // Next Case: we start again but before, we rank the points
    Compute      = %F;
  case 21 then
    // We replace W by C
    x_Worth      = x_Contr;
    f_Worth      = f_Contr;
    data_next(1) = f_Worth;
    data_next(4) = x_Worth;
    data_next(7) = 22; // Next Case: we start again but before, we rank the points
    Compute      = %F;
  case 22 then
    // Ranking the points
    f_nm = [];
    x_nm = list(); disp(n)
    for i=1:n+1
      f_nm(i)   = data_next(25)(i);
      x_nm($+1) = data_next(26)(i);
    end
    f_nm(1)   = f_Worth;
    f_nm($-1) = f_Good;
    f_nm($)   = f_Best;

    x_nm(1)   = x_Worth;
    x_nm($-1) = x_Good;
    x_nm($)   = x_Best;

    [f_sort, I] = sort(f_nm);

    f_Worth = f_nm(I(1));
    f_Good  = f_nm(I($-1));
    f_Best  = f_nm(I($));

    x_Worth = x_nm(I(1));
    x_Good  = x_nm(I($-1));
    x_Best  = x_nm(I($));

    // Updating the simplex
    for i=1:n+1
      data_next(25)(i) = f_nm(i);
      data_next(26)(i) = x_nm(i);
    end
    
    data_next(1) = f_Worth;
    data_next(2) = f_Good;
    data_next(3) = f_Best;
    data_next(4) = x_Worth;
    data_next(5) = x_Good;
    data_next(6) = x_Best;
    
    if (x_hist_is_defined) then
      data_next(23)($+1)  = list();
      data_next(23)($)(1) = x_Best;
      data_next(23)($)(2) = x_Good;
      data_next(23)($)(3) = x_Worth;
    end

    if (f_hist_is_defined) then
      data_next(22)($+1) = f_Best;
    end

    data_next(7) = 1;
    Compute      = %F;
  case 23 then
    // Computation of M
    x_next        = x_Median;
    data_next(21) = x_Median;
    data_next(7)  = 24;
    Compute       = %T;
  case 24 then
    // Computation of f(M)
    f_Median      = f_current;
    data_next(20) = f_Median;
    data_next(7)  = 20;
    eval_Func     = eval_Func + 1;
    data_next(24) = eval_Func;
    Compute       = %F;
  end
  
  data_current = data_next; // when in the while loop, we update the data_current
                            // structure to update the state of the Nelder and Mead
                            // method
end

if (f_hist_is_defined) then
  f_hist = data_next(22);
end
if (x_hist_is_defined) then
  x_hist = data_next(23);
end
endfunction


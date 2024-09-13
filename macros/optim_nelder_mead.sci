function [x_opt, x_history] = optim_nelder_mead(onm_f, x0, ItMX, Tol, MaxEvalFunc, kelley_restart, kelley_alpha)

[nargout, nargin] = argn();

if (size(x0,2)~=size(x0,1)+1) then
  printf('nelder_mead: you must give %d starting points\n',size(x0,1)+1);
end

n = size(x0,1);

x_history_defined = (nargout==2);

if (x_history_defined) then
  x_history = list();
  x_history($+1) = list();
  for i=1:n+1
    x_history($)($+1) = x0(:,i);
  end
end

if (~isdef('ItMX','local')) then
  ItMX = 100;
end
if (~isdef('Tol','local')) then
  Tol = 0.0;
end
if (~isdef('MaxEvalFunc','local')) then
  MaxEvalFunc = 10*ItMX;
end
if (~isdef('kelley_restart','local')) then
  kelley_restart = %F;
end
if (~isdef('kelley_alpha','local')) then
  kelley_alpha = 1e-4;
end

if ~isdef('onm_f','local') then
  error('optim_nelder_mead: onm_f is mandatory');
else
  if typeof(onm_f)=='list' then
    deff('y=_nm_f(x)','y=onm_f(1)(x, onm_f(2:$))');
  else
    deff('y=_nm_f(x)','y=onm_f(x)');
  end
end

// We compute and sort the objective function values of these 3 points

x_nm = list();

for i=1:n+1
  f_nm(i)   = _nm_f(x0(:,i));
  x_nm($+1) = x0(:,i);
end


[f_sort, I] = sort(f_nm);
f_Worth = f_nm(I(1));
f_Good  = f_nm(I($-1));
f_Best  = f_nm(I($));
x_Worth = x_nm(I(1));
x_Good  = x_nm(I($-1));
x_Best  = x_nm(I($));

if (Log) then
  printf('optim_nelder_mead: initialization\n');
  printf('f_Worth = %f f_Good = %f f_Best = %f\n', f_Worth, f_Good, f_Best);
end

FirstIteration = %T;
Iteration = 0;
f_Old     = f_Best;
eval_Func = 0;

// The main loop
//while(((eval_Func < MaxEvalFunc) & (Iteration<ItMX) & (abs((f_Best - f_Old)/f_Old)>Tol)) | (FirstIteration))
while(((eval_Func < MaxEvalFunc) & (Iteration<ItMX)) | (FirstIteration))
  f_Old          = f_Best;
  FirstIteration = %F;
  Iteration = Iteration + 1;
  
  x_Median = (x_Best + x_Good) / 2.0
  
  // First, we reflect the worth point
  if (Log) then
    printf('optim_nelder_mead: iteration = %d - reflexion\n', Iteration);
  end
  
  x_Reflexion = 2.0 * x_Median - x_Worth;  
  f_Reflexion = _nm_f(x_Reflexion);
  eval_Func   = eval_Func + 1;
    
  if (f_Reflexion<f_Good) then
    // Case I: Reflect or extend
    if (Log) then
      printf('optim_nelder_mead: reflect or extend - iteration = %d\n', Iteration);
    end
    if (f_Best<f_Reflexion) then
      f_Worth = f_Reflexion;
      x_Worth = x_Reflexion;
    else
      if (Log) then
        printf('optim_nelder_mead: extend\n');
      end
      x_Expansion = 2.0 * x_Reflexion - x_Median;
      f_Expansion = _nm_f(x_Expansion);
      eval_Func   = eval_Func + 1;
      if (f_Expansion<f_Best) then
        f_Worth = f_Expansion;
        x_Worth = x_Expansion;
      else
        f_Worth = f_Reflexion;
        x_Worth = x_Reflexion;
      end  
    end
  else
    // Case II: Extend or Shrink
    if (Log) then
      printf('optim_nelder_mead: extend or shrink - iteration = %d\n', Iteration);
    end

    if (f_Reflexion<f_Worth) then
      f_Worth = f_Reflexion;
      x_Worth = x_Reflexion;
    end
    if (Log) then
      printf('optim_nelder_mead: extend\n');
    end
    x_C1 = (x_Median + x_Worth) / 2.0;
    x_C2 = (x_Median + x_Reflexion) / 2.0;
    f_C1 = _nm_f(x_C1);
    eval_Func = eval_Func + 1;
    f_C2 = _nm_f(x_C2);
    eval_Func = eval_Func + 1;
    if (f_C1<f_C2) then
      f_Contr = f_C1;
      x_Contr = x_C1;
    else
      f_Contr = f_C2;
      x_Contr = x_C2;
    end
    if (f_Contr<f_Worth) then
      f_Worth = f_Contr;
      x_Worth = x_Contr;
    else
      if (Log) then
        printf('optim_nelder_mead: shrink - iteration = %d\n', Iteration);
      end
      x_Shrink  = (x_Best + x_Worth) / 2.0;
      f_Shrink  = _nm_f(x_Shrink);
      eval_Func = eval_Func + 1;
      f_Worth   = f_Shrink;
      x_Worth   = x_Shrink;
      
      f_Median  = _nm_f(x_Median);
      eval_Func = eval_Func + 1;
      f_Good    = f_Median;
      x_Good    = x_Median;
    end
  end

  // Ranking the points
  
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

  if (kelley_restart) then
    // Kelley restart
    // Computation of the simplex gradient
    df = f_nm(2:$) - f_nm(1);
    for i=2:length(x_nm)
      V  = x_nm(i) - x_nm(1);
      Df(i-1) = V'*df;
    end
    // Test
    if (f_Best - f_Old > -kelley_alpha * norm(Df)^2) then
      if (Log) then
        printf('optim_nelder_mead: Kelley restart\n');
      end
      // Simplex restart
      for i=1:n-1
        Aux(i) = norm(x_nm(1)-x_nm(i))
      end
      sigma_m_V = min(Aux);
      _beta = 0.5*sigma_m_V*sign(Df);
      y_nm = list();
      y_nm(1) = x_nm(1);
      for i=1:length(x_nm(1))
        y_nm(i+1)    = x_nm(1);
        y_nm(i+1)(i) = y_nm(i+1)(i) + _beta(i);
      end
      x_nm = list();
      for i=1:length(y_nm)
        x_nm(i) = y_nm(i);
        f_nm(i) = _nm_f(x_nm(i));
        eval_Func = eval_Func + 1;
      end
      y_nm = list();
    end
  end
  
  if (Log) then
    printf('optim_nelder_mead: f_Worth = %f f_Good = %f f_Best = %f\n', f_Worth, f_Good, f_Best);
  end

  if (x_history_defined) then
    x_history($+1) = list();
    for i=1:n+1
      x_history($)($+1) = x_nm(i);
    end
  end
end

x_opt = x_Best;
endfunction

function [x_opt, x_history] = optim_LFOP(lfop_df, x0, delta_t, delta_step, m, delta_inc, GradTOL, MaxEvalFunc, Log, p_start)
// df  : the gradient of the objective function. Can be a list like list(df, param1, param2)
// x0      : the initial solution
// delta_t : the temporal step size for the resolution of the dynamical system
// delta_step  : maximal allowable step size
// m           : id ak+1'.ak is non positive during m iterations then delta_t = delta_t / 2
//               and we restart from (xk+xk+1) / 2
// delta_inc   : increase factor of the dilatation coeff of delta_t
// GradTOL     : we stop the algorithm if the gradient is below the GradTOL level.
// MaxEvalFunc : we stop the algorithm when we have computed more than MaxEvalFunc the gradient.
// Log         : if %T, we display some informations during the running of the algorithm
// p_start     : starting value of the dilatation coefficient: 1.01 by default
// x_opt       : the last xk computed (the best solution)
// x_history   : the list of xk computed during the run of the optimization method

// Option for compilation of the deff's
option = 'n'; // 'c' for compiled deff. 'n' is better for debugging

[nargout, nargin] = argn();

if ~isdef('delta_step','local') then
  delta_step = 1;
end
if ~isdef('m','local') then
  m = 5;
end
if ~isdef('delta_inc','local') then
  delta_inc = 0.001;
end
if ~isdef('GradTOL','local') then
  GradTOL = 1e-5;
end
if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('p_start','local') then
  p_start  = 1.01;
end
if ~isdef('lfop_df','local') then
  error('optim_LFOP: lfop_df is mandatory');
else
  if typeof(lfop_df)=='list' then
    deff('y=_lfop_df(x)','y=lfop_df(1)(x, lfop_df(2:$))', opt=option);
  else
    deff('y=_lfop_df(x)','y=lfop_df(x)', opt=option);
  end
end
if ~isdef('x0','local') then
  error('optim_LFOP: x0 is mandatory');
end
if ~isdef('delta_t','local') then
  // We use a formula to size automatically the value of delta_t
  delta_t = sqrt(delta_step / (5*norm(_lfop_df(x0))));
end

Debug = %F;
Exit  = %F;
Cell     = 1;
EvalFunc = 0;
xk     = x0;
xk_m_1 = x0;
xk_p_1 = x0;
xk_p_2 = x0;

vk     = [];
vk_m_1 = [];
vk_p_1 = [];

a0 = [];
ak = [];
ak_p_1 = [];

x_history_defined = (nargout==2);
if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end
while((EvalFunc<MaxEvalFunc)&~Exit)
  select Cell
  case 1  then
    // Step 1 - Cell 1
    // First step initialisation
    i = 0;
    j = 2;
    s = 0;
    p = p_start; // 1 avant
    k = -1;
    Cell = 2;
  case 2  then
    // Step 2 - Cell 2
    // Second step initialisation
    a0 = - _lfop_df(x0);
    ak     = a0;
    ak_p_1 = a0;
    EvalFunc = EvalFunc + 1;
    v0 = 0.5*a0*delta_t;
    vk = v0;
    vk_m_1 = v0;
    vk_p_1 = v0;
    
    Cell = 3;
  case 3  then
    // Step 3 - Cell 3
    k = k + 1;
    
    xk_m_1 = xk;
    xk     = xk_p_1;
    vk_m_1 = vk;
    vk     = vk_p_1;
    ak     = ak_p_1;
    
    Cell = 4;
  case 4  then
    // Step 4 - Cell 4;
    if (Debug) then
      printf('||vk||delta_t = %f delta_t = %f p = %f delta_step = %f\n', norm(vk)*delta_t, delta_t, p, delta_step);
    end
    if (Log) then
      printf('EvalFunc = %d delta_t = %f p = %f norm(vk) = %f\n', EvalFunc, delta_t, p, norm(vk));
    end
    if (norm(vk)*delta_t < delta_step) then
      Cell = 5;
    else
      vk = delta_step*vk / (delta_t*norm(vk));
      Cell = 6;
    end
  case 5  then
    // Step 5b - Cell 5
    p = p + delta_inc;
    delta_t = p * delta_t;
    Cell = 6;
  case 6  then
    // Step 5a - Cell 6
    if (s<m) then
      Cell = 8;
    else
      Cell = 7;
    end
  case 7  then
    // Step 5b - Cell 7
    delta_t = delta_t / 2;
    xk = (xk + xk_m_1) / 2;
    vk = (vk + vk_m_1) / 4;
    s  = 0;
    Cell = 8;
  case 8  then
    // Step 5 - Cell 8
    xk_p_1 = xk + vk*delta_t;
    Cell = 9;
  case 9  then
    // Step 6 - Cell 9
    ak_p_1 = - _lfop_df(xk_p_1);
    
    if (x_history_defined&~isempty(xk_p_1)) then
      x_history($+1) = xk_p_1;
    end
    
    EvalFunc = EvalFunc + 1;
    vk_p_1 = vk + ak_p_1*delta_t;
    Cell = 10;
  case 10 then
    // Step 7a - Cell 10
    if (Debug) then
      printf('ak_p_1''*ak ='); disp(ak_p_1'*ak);
    end
    if (ak_p_1'*ak > 0) then
      Cell = 11;
    else
      Cell = 12;
    end
  case 11 then
    // Step 7a - Cell 11
    s = 0;
    Cell = 13;
  case 12 then
    // Step 7a - Cell 12
    s = s + 1;
    p = p_start;
    Cell = 13;
  case 13 then
    // Step 7 - Cell 13
    if (norm(ak_p_1)<=GradTOL) then
      Cell = 14;
    else
      Cell = 15;
    end
  case 14 then
    // Stop Cell;
    Exit = %T;
  case 15 then
    // Step 8 - Cell 15
    if (norm(vk_p_1)>norm(vk)) then
      Cell = 16;
    else
      Cell = 17;
    end
  case 16 then
    // Step 8 - Cell 16
    i = 0;
    Cell = 3; // Go to step 3
  case 17 then
    // Step 8 - Cell 17
    xk_p_2 = (xk_p_1 + xk) / 2;
    i = i + 1;
    Cell = 18;
  case 18 then
    // Step 9 - Cell 18
    if (i<=j) then
      Cell = 19;
    else
      Cell = 20;
    end
  case 19 then
    // Step 9 - Cell 19
    vk_p_1 = (vk_p_1 + vk)/4;
    k = k + 1;
    xk_m_1 = xk;
    xk     = xk_p_1;
    xk_p_1 = xk_p_2;
    vk     = vk_p_1;
    vk_m_1 = vk;
    ak     = ak_p_1;
    
    Cell = 9; // Go to step 6
  case 20 then  
    // Step 9 - Cell 20
    vk_p_1 = zeros(size(x0,1),size(x0,2));
    j = 1;
    k = k + 1;
    xk_m_1 = xk;
    xk     = xk_p_1;
    xk_p_1 = xk_p_2;
    vk     = vk_p_1;
    vk_m_1 = vk;
    ak     = ak_p_1;
    
    Cell = 9; // Go to step 6
  end
end
x_opt = xk;
endfunction

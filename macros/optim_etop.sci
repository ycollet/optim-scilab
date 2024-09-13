function [x_opt, x_history] = optim_etop(etop_df, x0, dt, ItMX, cgmethod, deltax_max, XTol, GTol, Log, Restart)
[nargout, nargin] = argn();

x_history_defined = (nargout==2);
if x_history_defined then
  x_history = list();
  x_history($+1) = x0;
end

if ~isdef('dt','local') then
  dt = 0.5;
end
if ~isdef('ItMX','local') then
  ItMX = 1000;
end
if ~isdef('Log','local') then
  Log = %F;
end
if ~isdef('XTol','local') then
  XTol = 0;
end
if ~isdef('GTol','local') then
  GTol = 0;
end
if ~isdef('cgmethod','local') then
  cgmethod = 'fr'; // pr for the other choice
end
if ~isdef('Restart','local') then
  Restart = length(x0)+1;
end
if ~isdef('deltax_max','local') then
  deltax_max = ones(size(x0,1),size(x0,2))*%inf;
end

deltax_max = abs(deltax_max);

xk   = x0;
xk_1 = xk;
_beta = 0;

dk   = etop_df(xk);
dk_1 = dk;
vk   = -dk*dt;
vk_1 = -dk*dt;
dt_1 = dt;

for i=1:ItMX
  if Log then
    printf('optim_etop: Iteration %d / %d \n', i, ItMX);
  end
  
  xk_1   = xk;
  dk_1   = dk;
  vk_1   = vk;
  
  xk = xk_1 + vk*dt^2;

 // if or(abs(xk-xk_1)>deltax_max) then
//    xk = xk_1 + deltax_max.*sign(xk - xk_1).*(abs(xk-xk_1)>deltax_max) + (xk - xk_1).*(~(abs(xk-xk_1)>deltax_max));
//  end

  dk = etop_df(xk);
  
  _deltaF = -(vk*dt)'*(-dk-dk_1)*dt*0.5; // Formule pas identique dans le programme ... a revoir
  _rho    = (xk - xk_1)'*dk_1;           // Formule identique dans le programme
  _theta  = - _rho / (_deltaF - _rho);   // Formule identique dans le programme

  // Heuristics for _theta - First part
  if (abs(_rho-_deltaF)<=100*%eps) then 
    _theta = 1; 
    if Log then
      printf('optim_etop: theta set to 1\n');
    end
  elseif _deltaF < 0 & _theta < 0 then
    _theta = - _theta;
    if Log then
      printf('optim_etop: mirroring theta\n');
    end
  end

  //_theta = min([1 max([_theta -1])]);
  
  xk_star = xk_1 + _theta*(xk - xk_1)*0.5;
  
//  if or(abs(xk_star-xk_1)>deltax_max) then
//    xk_star = xk_1 + deltax_max.*sign(xk_star - xk_1).*(abs(xk_star-xk_1)>deltax_max) + (xk_star - xk_1).*(~(abs(xk_star  -xk_1)>deltax_max));
//  end

  xk = xk_star;
  
  dk = etop_df(xk);

  if cgmethod=='fr' then
    _beta = norm(dk)^2 / norm(dk_1)^2;
  elseif cgmethod=='pr' then
    _beta = (dk - dk_1)'*dk / norm(dk_1)^2;
  else
    error('optim_etop: cgmethod must be ''pr'' or ''fr''\n');
  end

  if modulo(i,Restart)==0 then
    _beta = 0;
    if Log then
      printf('optim_etop: conjugate gradient restart\n');
    end
  end
  
  _beta = min([_beta 1]);
  vk    = (-dk + _beta*vk_1/dt_1)*dt;
  
  if x_history_defined then
    x_history($+1) = xk;
  end

  if Log then
    printf('optim_etop: dt = %f, ||grad f|| = %f ||vk|| = %f dF = %f theta = %f rho = %f beta = %f\n', dt, norm(dk), norm(vk), _deltaF, _theta, _rho, _beta);
    printf('optim_etop: norm(xk-xk_1) = %f norm(dk) = %f\n', norm(xk-xk_1), norm(dk));
  end

  dt_1 = dt;
  
  // Heuristics for _theta - Second part
  if _deltaF > 0 & _theta > 0 then
    if Log then
      printf('optim_etop: contraction\n');
    end
    dt = dt * 0.5;
  elseif _deltaF <= 0 & _theta >= 0 then
    if Log then
      printf('optim_etop: dilatation\n');
    end
    dt = dt * 1.5;
  elseif _deltaF < 0 then
    if Log then
      printf('optim_etop: dilatation\n');
    end
    dt = dt * 1.5;
  end

  // Convergence tests
  if norm(xk - xk_1)<XTol then
    if Log then
      printf('optim_etop: break on XTol\n');
    end
    break;
  end
  if norm(dk)<GTol then
    if Log then
      printf('optim_etop: break on GTol\n');
    end
    break;
  end
  if i==ItMX then
    if Log then
      printf('optim_etop: break on iteration reached\n');
    end
  end
end
x_opt = xk;
endfunction

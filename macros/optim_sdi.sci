function [x_opt, x_history, x_pool_history] = optim_sdi(sdi_f, x0, noise_type, noise_param, target, weight, NbLoops, NbPoints)

[nargout, nargin] = argn();

x_history_defined = (nargout>=2);
x_pool_history_defined = (nargout==3);

if ~isdef('sdi_f','local') then
  error('optim_sdi: sdi_f is mandatory');
else
  if typeof(sdi_f)=='list' then
    deff('y=_sdi_f(x)','y=sdi_f(1)(x, sdi_f(2:$))');
  else
    deff('y=_sdi_f(x)','y=sdi_f(x)');
  end
end
if ~isdef('x0','local') then
  error('optim_sdi: x0 is mandatory');
end
if ~isdef('noise_type','local') then
  for i=1:length(x0)
    noise_type(i) = 'uniform';
  end
end
if ~isdef('noise_param','local') then
  noise_param = [ones(length(x0),1); zeros(length(x0),1)];
end
f_value = _sdi_f(x0);
if ~isdef('target','local') then
  target = zeros(size(f_value,1), size(f_value,2));
end
if ~isdef('weight','local') then
  weight = ones(size(f_value,1), size(f_value,2));
end
if ~isdef('NbLoops','local') then
  NbLoops = 10;
end
if ~isdef('NbPoints','local') then
  NbPoints = 10;
end

if (x_history_defined) then
  x_history = list();
  x_history($+1) = x0;
end
if (x_pool_history_defined) then
  x_pool_history = list();
end

for i=1:NbLoops
  // Generating a pool of points
  for j=1:NbPoints
    for k=1:length(x0)
      XPool(k,j) = noise_param(1,k)*rand(1,1,noise_type(k)) + noise_param(2,k) + x0(k);
    end
  end
  if (x_pool_history_defined) then
    x_pool_history($+1) = list();
    for j=1:NbPoints
      x_pool_history($)($+1) = XPool(:,j);
    end
  end
  // Computing the objective functions values
  for j=1:NbPoints
    f_value = _sdi_f(XPool(:,j));
    for k=1:length(f_value)
      FPool(k,j) = f_value(k);
    end
  end
  // Computing the distance to the target
  for j=1:NbPoints
    DistPool(j) = sqrt(sum(weight.*(FPool(:,j) - target).^2));
  end
  // Finding the best point
  [Value, Index] = min(DistPool);
  // Updating x0
  x0 = XPool(:,Index(1));
  if (x_history_defined) then
    x_history($+1) = x0;
  end
end
x_opt = x0;
endfunction

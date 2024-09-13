function plot_d(f, x0, d, lbounds, ubounds, restrtomin, hlow, hup)
// function plot_d: this function plot a function of several variables in the direction given by vector d
//                  between lbounds and ubounds
// f          : the function to plot
// x0         : the initial point
// d          : the direction
// lbounds    : a vector of lower boundary for the definition domain of the function f
// ubounds    : a vector of upper boundary for the definition domain of the function f
// restrtomin : a boolean to plot only up to the min of the function (optional parameter)
// hlow       : lower boundary on the h parameter (optional parameter)
// hup        : upper boundary on the h parameter (optional parameter)

[nargout, nargin] = argn();

if ((nargin==5)|(nargin==6)) then 
  HLowBounded = %F;
  HUpBounded  = %F;
elseif (nargin==7) then
  HLowBounded = %T;
  HUpBounded  = %F;
elseif (nargin==8) then
  HLowBounded = %T;
  HUpBounded  = %T;
else
  error('plot_d: wrong number of parameters');
  return;
end

// Normalization of the search direction

I = find(abs(d)<=%eps);
d(I) = %eps;

d_norm = d/norm(d);

// Searching for h_min, h_max

h_max_vec = (ubounds - x0)./abs(d_norm);
h_min_vec = (x0 - lbounds)./abs(d_norm);

h_min = min(h_min_vec); h_min = - h_min;
h_max = min(h_max_vec);

if (HLowBounded) then h_min = hlow; end
if (HUpBounded)  then h_max = hup;  end

// Plot the function

h_plot_d = h_min:(h_max-h_min)/100.0:h_max;
for i=1:size(h_plot_d,2)
  y_plot_d(i) = f(x0+h_plot_d(i)*d_norm);
end
[m, I] = min(y_plot_d);
if (I(1)<size(y_plot_d,2)) then I(1) = I(1) + 1; end

plot(h_plot_d(1:I(1)), y_plot_d(1:I(1)), 'k-');
xtitle('d-cutting plane of function f','h','f(x0+h*d)');
endfunction

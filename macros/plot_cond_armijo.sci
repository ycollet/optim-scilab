function plot_cond_armijo(xmin, xstep, xmax, f, fs, c1)
// Function to plot the armijo condition
// xmin  : starting point
// xstep : plot every xstep points
// xmax  : final point
// f     : objective function (must be 1D)
// fs    : derivative function (must be 1D)
// c1    : parameter of the armijo condition

color_1 = color('red');
color_2 = color('green');

color_list_1 = [];
color_list_2 = [];
y = [];
ys = [];
ys_wc = [];

fk  = f(xmin);
fsk = fs(xmin);
x = xmin:xstep:xmax;
y_arm = fk + c1*(x-xmin)*fsk;

for i=1:size(x,2)
  y(i)     = f(x(i));
  ys(i)    = fs(x(i));
  
  if (y(i) <= y_arm(i)) then
    color_list_1(i) = color_1;
  else
    color_list_1(i) = color_2;
  end
end

I1 = find(color_list_1==color_2);

yaux         = y;
yaux_arm     = y_arm;
yaux(I1)     = %nan;
yaux_arm(I1) = %nan;

plot2d(x, yaux,     color_1);
plot2d(x, yaux_arm, color_1);

I2 = find(color_list_1==color_1);
I2 = I2(1,1:$-1);

yaux     = y;
yaux_arm = y_arm;

yaux(I2)     = %nan;
yaux_arm(I2) = %nan;


plot2d(x, yaux,     color_2);
plot2d(x, yaux_arm, color_2);

xtitle('Bound Armijo condition','x','f(x)');
legends(['Condition verified','Condition not verified'],[color_1 color_2;1 1], opt='ur');
endfunction

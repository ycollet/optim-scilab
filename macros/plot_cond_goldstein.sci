function plot_cond_goldstein(xmin, xstep, xmax, f, fs, c)
// Function to plot the goldstein condition
// xmin  : starting point
// xstep : plot every xstep points
// xmax  : final point
// f     : objective function (must be 1D)
// fs    : derivative function (must be 1D)
// c1    : fist parameter of the wolf condition
// c2    : second parameter of the wolf condition

color_1 = color('red');
color_2 = color('green');

color_list_1 = [];
color_list_2 = [];
y = [];
ys = [];
ys_wc = [];

fk    = f(xmin);
fsk   = fs(xmin);
x     = xmin:xstep:xmax;
y1_wc = fk + (1-c)*(x-xmin)*fsk;
y2_wc = fk + c*(x-xmin)*fsk;

for i=1:size(x,2)
  y(i)  = f(x(i));
  ys(i) = fs(x(i));
  
  if (y(i) >= y1_wc(i)) then
    color_list_1(i) = color_1;
  else
    color_list_1(i) = color_2;
  end

  if (y(i)<=y2_wc(i)) then
    color_list_2(i) = color_1;
  else
    color_list_2(i) = color_2;
  end
end

subplot(2,1,1);
I1 = find(color_list_1==color_2);
yaux        = y;
yaux_wc     = y1_wc;

yaux(I1)    = %nan;
yaux_wc(I1) = %nan;

plot2d(x, yaux,    color_1);
plot2d(x, yaux_wc, color_1);

I2 = find(color_list_1==color_1);
I2 = I2(1,1:$-1);

yaux        = y;
yaux_wc     = y1_wc;

yaux(I2)    = %nan;
yaux_wc(I2) = %nan;

plot2d(x, yaux,    color_2);
plot2d(x, yaux_wc, color_2);

xtitle('Bound Goldstein condition','x','f(x)');
legends(['Condition verified','Condition not verified'],[color_1 color_2;1 1], opt='ur');

subplot(2,1,2);

I1 = find(color_list_2==color_2);
yaux        = y;
yaux_wc     = y2_wc;
yaux(I1)    = %nan;
yaux_wc(I1) = %nan;

plot2d(x, yaux,    color_1);
plot2d(x, yaux_wc, color_1);

I2 = find(color_list_2==color_1);
I2 = I2(1,2:$);

yaux        = y;
yaux_wc     = y2_wc;
yaux(I2)    = %nan;
yaux_wc(I2) = %nan;

plot2d(x, yaux,    color_2);
plot2d(x, yaux_wc, color_2);

xtitle('Slope Goldstein condition','x','f(x)');
legends(['Condition verified','Condition not verified'],[color_1 color_2;1 1],opt='ur');
endfunction

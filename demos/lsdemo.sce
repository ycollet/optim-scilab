Log      = 1;
TOL      = 1e-4;
ItMAX    = 50;
h        = 4.0;
x0       = 0;
// testfunc = 'more_thuente_1';
// testfunc = 'more_thuente_2';
testfunc = 'more_thuente_3'; // Pas bon pour newton ni pour secant
// testfunc = 'YOK_1';
// testfunc = 'YOK_2';
// testfunc = 'YOK_3'; // Pas bon pour newton

// testmethod = 'dicho';
// testmethod = 'goldsect';
testmethod = 'newton';
// testmethod = 'secant';
// testmethod = 'backtrack';
// testmethod = 'polynom';
// testmethod = '2pts_bracket';
// testmethod = '3pts_bracket';

//
// Parametrisation of the various test functions
//

if (testfunc=='more_thuente_1') then
  deff('y=df_more_thuente_1(x)','y=derivative(more_thuente_1,x)');
  deff('y=dff_more_thuente_1(x)','[tmp,y]=derivative(more_thuente_1,x)');
  f     = more_thuente_1;
  fs    = df_more_thuente_1;
  fss   = dff_more_thuente_1;
  xmin  = 0;
  xmax  = 16;
  Title = 'More et Thuente 1';
end
if (testfunc=='more_thuente_2') then
  deff('y=df_more_thuente_2(x)','y=derivative(more_thuente_2,x)');
  deff('y=dff_more_thuente_2(x)','[tmp,y]=derivative(more_thuente_2,x)');
  f     = more_thuente_2;
  fs    = df_more_thuente_2;
  fss   = dff_more_thuente_2;
  xmin  = 0;
  xmax  = 2;
  Title = 'More et Thuente 2';
end
if (testfunc=='more_thuente_3') then
  deff('y=df_more_thuente_3(x)','y=derivative(more_thuente_3,x)');
  deff('y=dff_more_thuente_3(x)','[tmp,y]=derivative(more_thuente_3,x)');
  f     = more_thuente_3;
  fs    = df_more_thuente_3;
  fss   = dff_more_thuente_3;
  xmin  = 0;
  xmax  = 2;
  Title = 'More et Thuente 3';
end
if (testfunc=='YOK_1') then
  deff('y=df_YOK_1(x)','y=derivative(YOK_1,x)');
  deff('y=dff_YOK_1(x)','[tmp,y]=derivative(YOK_1,x)');
  f     = YOK_1;
  fs    = df_YOK_1;
  fss   = dff_YOK_1;
  xmin  = 0;
  xmax  = 2;
  Title = 'Yanai, Ozawa et Kaneko 1';
end
if (testfunc=='YOK_2') then
  deff('y=df_YOK_2(x)','y=derivative(YOK_2,x)');
  deff('y=dff_YOK_2(x)','[tmp,y]=derivative(YOK_2,x)');
  f     = YOK_2;
  fs    = df_YOK_2;
  fss   = dff_YOK_2;
  xmin  = 0;
  xmax  = 2;
  Title = 'Yanai, Ozawa et Kaneko 2';
end
if (testfunc=='YOK_3') then
  deff('y=df_YOK_3(x)','y=derivative(YOK_3,x)');
  deff('y=dff_YOK_3(x)','[tmp,y]=derivative(YOK_3,x)');
  f     = YOK_3;
  fs    = df_YOK_3;
  fss   = dff_YOK_3;
  xmin  = 0;
  xmax  = 2;
  Title = 'Yanai, Ozawa et Kaneko 3';
end

// We plot the selected test function

x = []; y = [];

x=xmin:(xmax-xmin)/100.0:xmax;

for i=1:size(x,2)
  y(i)=f(x(i));
end

plot(x,y,'b-');
xgrid(2);

if (testmethod=='newton') then
  printf('ls_newton:\n\n');
  [h_fin, history] = ls_newton(f, fs, fss, x0, -fs(x0), h, Log, TOL, ItMAX);
end

if (testmethod=='secant') then
  printf('ls_secant:\n\n');
  [h_fin, history] = ls_secant(f, fs, fss, x0, -fs(x0), h, Log, TOL, ItMAX);
end

if (testmethod=='dicho') then
  // The dichotomy process is very sensitive to the size step
  printf('ls_dicho:\n\n');
  [h_fin, history] = ls_dicho(f, fs, fss, x0, -fs(x0), h/1000, Log, TOL, ItMAX);
end

if (testmethod=='backtrack') then
  printf('ls_backtrack:\n\n');
  [h_fin, history] = ls_backtrack(f, fs, fss, x0, -fs(x0), h, Log, TOL, ItMAX);
end

if (testmethod=='polynom') then
  printf('ls_polynom:\n\n');
  [h_fin, history] = ls_polynom(f, fs, fss, x0, -fs(x0), h/1000, Log, TOL, ItMAX);
end

if (testmethod=='goldsect') then
  printf('ls_goldsect:\n\n');
  [h_fin, history] = ls_goldsect(f, fs, fss, x0, -fs(x0), h, Log, TOL, ItMAX);
end

if (testmethod=='2pts_bracket') then
  printf('ls_2pts_bracket:\n\n');
  [h_fin, history] = ls_2pts_bracket(f, fs, fss, x0, -fs(x0), h, Log, TOL, ItMAX);
  printf('fs(h_min) = %f - fs(h_max) = %f\n', fs(x0-h_fin(1)*fs(x0)), fs(x0-h_fin(2)*fs(x0)));
  printf('f(h_min)  = %f - f(h_max)  = %f\n', f(x0-h_fin(1)*fs(x0)), f(x0-h_fin(2)*fs(x0)));
  printf('h_min     = %f - h_max     = %f\n', h_fin(1), h_fin(2));
end

if (testmethod=='3pts_bracket') then
  printf('ls_3pts_bracket:\n\n');
  [h_fin, history] = ls_3pts_bracket(f, fs, fss, x0, -fs(x0), h, Log, TOL, ItMAX);
  printf('fs(h_min) = %f - fs(h_mid) = %f - fs(h_max) = %f\n', fs(x0-h_fin(1)*fs(x0)), fs(x0-h_fin(2)*fs(x0)), fs(x0-h_fin(3)*fs(x0)));
  printf('f(h_min)  = %f - f(h_mid)  = %f - f(h_max)  = %f\n', f(x0-h_fin(1)*fs(x0)), f(x0-h_fin(2)*fs(x0)), f(x0-h_fin(3)*fs(x0)));
  printf('h_min     = %f - h_mid     = %f - h_max     = %f\n', h_fin(1), h_fin(2), h_fin(3));
end

y = [];
x = [];

for i=1:length(history)
  y(i)   = f(history(i));
  x(:,i) = history(i);
end

xtitle(testmethod + ' - ' + testfunc,'alpha','f(x0+alpha.d)');
plot(x,y,'r+');



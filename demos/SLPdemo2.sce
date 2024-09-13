
delta_ml = 0.1;
upper    = ones(10,1);
lower    = zeros(10,1);
x0       = sort(rand(10,1));
ItMX        = 1000;
ItMX_SIP    = 10;
MaxEvalFunc = 1000;

// vv_open_obj
// vv_open_constr
deff('y=df_vv_open_obj(x)','y=derivative(vv_open_obj,x)''');
deff('y=df_vv_open_constr(x)','y=derivative(vv_open_constr,x)''');
[x_opt, x_history, ml_history] = optim_slp(vv_open_obj, ...
                                           df_vv_open_obj, ...
                                           vv_open_constr, ...
                                           df_vv_open_constr, ...
                                           x0, ItMX, MaxEvalFunc, delta_ml, upper, lower, ItMX_SIP);

disp(x_opt)

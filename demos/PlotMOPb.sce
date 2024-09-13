// binh_1         - not constrained - 2D
// fonseca_1      - not constrained - 2D
// fonseca_2      - not constrained - 2D
// kursawe_1      - not constrained - 2D
// laumanns_1     - not constrained - 2D
// lis_1          - not constrained - 2D
// murata_1       - not constrained - 2D
// poloni_1       - not constrained - 2D
// quagliar_1     - not constrained - 2D
// rendon_1       - not constrained - 2D
// rendon_2       - not constrained - 2D
// schaffer_1     - not constrained - 2D
// schaffer_2     - not constrained - 2D
// binh_2         - not constrained - 3D
// viennet_1      - not constrained - 3D
// viennet_2      - not constrained - 3D
// viennet_3      - not constrained - 3D
// deb_1          - not constrained - 2D
// deb_2          - not constrained - 2D
// deb_3          - not constrained - 2D
// deb_4          - not constrained - 2D
// deb_5          - not constrained - 2D
// deb_6          - not constrained - 2D
// deb_7          - not constrained - 2D
// deb_8          - not constrained - 2D
// deb_8_bis      - not constrained - 2D
// veldmop_1      - not constrained - 2D
// veldmop_2      - not constrained - 2D
// veldmop_3      - not constrained - 2D
// veldmop_4      - not constrained - 2D
// veldmop_5      - not constrained - 3D
// veldmop_6      - not constrained - 2D
// veldmop_7      - not constrained - 3D
// meca_1         - not constrained - 2D
// trigo_1        - not constrained - 2D
// trigo_2        - not constrained - 2D
// trigo_3        - not constrained - 2D
// trigo_4        - not constrained - 2D
// trigo_5        - not constrained - 2D
// trigo_5_bis    - not constrained - 2D
// trigo_6        - not constrained - 2D
// trigo_6_bis    - not constrained - 2D
// trigo_7        - not constrained - 2D
// trigo_8        - not constrained - 2D
// trigo_9        - not constrained - 2D

// belegundu_1    - constrained     - 2D
// binh_3         - constrained     - 2D
// jimenez_1      - constrained     - 2D
// kita_1         - constrained     - 2D
// obayashi_1     - constrained     - 2D
// osyczka_1      - constrained     - 2D
// osyczka_2      - constrained     - 2D
// srinivas_1     - constrained     - 2D
// tanaka_1       - constrained     - 2D
// binh_4         - constrained     - 3D
// tamaki_1       - constrained     - 3D
// viennet_4      - constrained     - 3D
// coello_1       - constrained     - 2D
// coello_3       - constrained     - 2D
// veldmopc_1     - constrained     - 2D
// veldmopc_2     - constrained     - 2D
// veldmopc_3     - constrained     - 3D
// hanne_1        - constrained     - 2D
// hanne_2        - constrained     - 2D
// hanne_3        - constrained     - 2D
// hanne_4        - constrained     - 2D
// hanne_5        - constrained     - 2D
// meca_2         - constrained     - 2D

//funcname = 'quagliar_1';
//funcname = 'rendon_1';
//funcname = 'rendon_2';
//funcname = 'schaffer_1';
//funcname = 'schaffer_2';
//funcname = 'binh_2';
//funcname = 'viennet_1';
//funcname = 'viennet_2';
//funcname = 'viennet_3';
//funcname = 'deb_1';
//funcname = 'deb_2'; 
//funcname = 'deb_3';
//funcname = 'deb_4';
//funcname = 'deb_5'; 
//funcname = 'deb_6';
funcname = 'deb_7';
//funcname = 'deb_8';
//funcname = 'veldmop_3';
//funcname = 'veldmop_4';
//funcname = 'veldmop_5';
//funcname = 'veldmop_6';
//funcname = 'veldmop_7';
//funcname = 'meca_1';
//funcname = 'trigo_1';
//funcname = 'trigo_2';
//funcname = 'trigo_3';
//funcname = 'trigo_4';
//funcname = 'trigo_5';
//funcname = 'trigo_5_bis';
//funcname = 'trigo_6';
//funcname = 'trigo_6_bis';
//funcname = 'trigo_7';
//funcname = 'trigo_8';
//funcname = 'trigo_9';
//funcname = 'belegundu_1';
//funcname = 'binh_3';
//funcname = 'jimenez_1';
//funcname = 'kita_1';
//funcname = 'osyczka_1';
//funcname = 'osyczka_2';
//funcname = 'srinivas_1';
//funcname = 'tanaka_1';
//funcname = 'binh_4';
//funcname = 'viennet_4';
//funcname = 'coello_1';
//funcname = 'coello_3'; // ?? Y contient %nan !!
//funcname = 'veldmopc_1';
//funcname = 'veldmopc_2';
//funcname = 'hanne_1';
//funcname = 'hanne_2';
//funcname = 'hanne_3';
//funcname = 'hanne_4';
//funcname = 'hanne_5';
//funcname = 'meca_2';

N_pts = 40000;

Min = eval('get_min_bound_'+funcname+'()');
Max = eval('get_max_bound_'+funcname+'()');
Aux = eval(funcname+'(Min)');

Y = zeros(N_pts,length(Aux));
tmp = eval(funcname+'(Min)');
Is2D = (length(tmp)==2);
Constrained = isdef('constr_'+funcname);

wID = waitbar('plotting points');
for i=1:N_pts
  if (modulo(i,100)==0) then
    waitbar(i/N_pts,wID);
  end
  Aux = (Max - Min).*rand(size(Min,1),size(Min,2)) + Min;
  if Constrained then
    while(or(eval('constr_'+funcname+'(Aux)')>0))
      Aux = (Max - Min).*rand(size(Min,1),size(Min,2)) + Min;
    end
  end
  Y(i,:) = eval(funcname+'(Aux)');
end
winclose(wID);

if Constrained then
  Title = sprintf('%s function - constrained', funcname);
else
  Title = sprintf('%s function - not constrained', funcname);
end

drawlater;
if Is2D then
  plot(Y(:,1),Y(:,2),'k.');
  xtitle(Title,'f1','f2');
  h = gca();
  h.children.children.mark_size = 1;
else
  param3d1(Y(:,1),Y(:,2),list(Y(:,3),0));
  xtitle(Title,'f1','f2','f3');
  h = gca();
  h.children.mark_size = 0.5;
  h.cube_scaling       = 'on';
end
drawnow;


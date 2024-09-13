deff('y=f(x)','y=x^2');
deff('y=df(x)','y=2*x');

x = -2:0.01:2;
h = log(1e-6):abs(log(1)-log(1e-6))/10:log(1);

ResFD = [];
for i=1:length(h)
  ResFD(i) = -%inf;
  for j=1:length(x)
    Aux = df(x(j)) - (f(x(j) + exp(h(i))) - f(x(j))) / exp(h(i));
    if (Aux > ResFD(i)) then ResFD(i) = Aux; end
  end
end

ResCFD = [];
for i=1:length(h)
  ResCFD(i) = -%inf;
  for j=1:length(x)
    Aux = df(x(j)) - (f(x(j) + exp(h(i))) - f(x(j) - exp(h(i)))) / (2*exp(h(i)));
    if (Aux > ResCFD(i)) then ResCFD(i) = Aux; end
  end
end

plot(h, ResFD, 'r-');
plot(h, ResCFD, 'g-');
xtitle('Comparaison dérivée / dérivée estimée','h','erreur');
legend(['FD','Central FD']);


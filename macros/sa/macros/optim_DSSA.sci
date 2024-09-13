function [xmin,x_history] = optim_DSSA(f,U,L,ItMX)
// *******************************************************
// optim_DSSA = Direct Search Simulated Annealing 
// The DSSA Code is created on Oct. 16, 2001 by A. Hedar 
// Modified on Aug. 13, 2005
// Reference: 
//       A. Hedar and M. Fukushima,Hybrid simulated annealing and direct   
//       search method for nonlinear unconstrained global optimization, 
//       Optimization Methods and Software, 17 (2002) 891-912.
// *******************************************************
// Inputs:
//  n    = Number of the objective variables   
//  U    = n-Dimension Vector of upper limits of the objective variables
//  L    = n-Dimension Vector of lower limits of the objective variables
//  f    = Objective function
//  ItMX = maximum number of iterations
//
// Outputs:
//  xmin   = Best point obtained by SAHPS 
//  fmin   = Best function value obtained by SAHPS
//  fcount = Number of function evaluations                   
// 

[nargout, nargin] = argn();

//
//Initialize counters
//

fcount = 0;	//number of functional evaluations
itc    = 0;	//number of iterations
x_history_defined = (nargout==2);

//
//Set the parameters	
//

//
// Generation of the initial point
//
rang    = max(U-L);
_center = (U+L)/2;
n       = length(U);
x       = zeros(n,1);
for i=1:n;
  x0(i)=L(i)+rand()*(U(i)-L(i));
end
x1 = 2*_center-x0;

//
//Set the parameters	
//

if nargin < 5 then ItMX=50*n; end

rd     = 0.5;    //reduction factor of temp
maxrit = n;      //number of iterations before reducte the temp
tol    = 1.0D-8; //termanation acurcy
nbest  = n;      //number of best stored functional values
rad    = 0.001;  //the redius of the balls in the best points List
mloopc = 2;		    //the number of iterating the main loop
edger  = 0.5;	   //reduction factor of edge length for simplex

//
//creating the initial simplex
//

x      = zeros(n,n+1);
x(:,1) = x0;
id     = eye(n,n);
edgel  = min(0.1*max(U-L),4);
edgebest = edgel;
df     = 0;
fv(1)  = f(x(:,1));
fcount = fcount+1;
edge   = edgel*ones(n,1);
for j=2:n+1; x(:,j)=x0+edge(j-1)*id(:,j-1); end
for j=2:n+1; fv(j)=f(x(:,j)); end;
fcount = fcount+n;
for j=1:n; d(j)=abs(fv(1)-fv(j+1)); end
df = norm(d);

if x_history_defined then
  x_history = list();
  x_history($+1) = list();
  for ii=1:size(x,2)
    x_history($)($+1) = x(:,ii);
  end
end

while(df<10 & edgel<4)
  edgel = 2*edgel;
  edge  = edgel*ones(n,1);
  for j=2:n+1; x(:,j) = x0+edge(j-1)*id(:,j-1); end
  for j=2:n+1; fv(j)  = f(x(:,j)); end;
  fcount = fcount+n;
  for j=1:n; d(j)=abs(fv(1)-fv(j+1)); end
  df = norm(d);          
end

//
//Order the vertices for the first time
//

xtmp    = zeros(n,n+1); 
[fs,is] = sort(-fv); fs = -fs; xtmp=x(:,is); x=xtmp; fv=fs;
maxdf   = fv(n+1)-fv(1);
itc     = itc+1;
// Tmax=(fv(n+1)-fv(1))/log(1/0.9);
// Tmin=0.00001*Tmax;
Tmax = max(0.001,(fv(n+1)-fv(1)))/log(1/0.9);
Tmin = 0.00001*Tmax;
T    = Tmax;
for j=1:nbest; bestxl(:,j)=x(:,1); bestfl(j)=fv(1);end

//
// Main loop
//

for loopc=1:mloopc;
  while(itc < ItMX & T > Tmin)  
    rditc = 0;
    while(rditc < maxrit)
      r      = 1; // number of the reflecting vertices
      accept = 0;
      while(accept == 0)
        p     = n-r+1;
        one2p = 1:p;
        xbar  = sum(x(:,one2p),2)/p; //centroid of p best point
        for j=1:r
          k = j+p;
          xr(:,j) = xbar+(0.2*rand()+0.9)*(xbar-x(:,k));
        end
        for j=1:r; fr(j) = f(xr(:,j)); end
        fcount    = fcount+r;
        [frs,irs] = sort(-fr); frs = -frs; ii=irs(1);
        frmin = frs(1); xrmin=xr(:,ii);
        prob  = exp(-(frmin-fv(1))/T);
        if frmin < fv(1) then
          for j=1:r; k=j+p; x(:,k)=xr(j); end
          bestla = 0;
          jj     = 1;
          while(bestla~=1 & jj < nbest)
            v(:,jj) = xrmin-bestxl(:,jj);
            if norm(v(:,jj)) < rad then
              bestla = 1;
            end
            jj = jj+1;
          end
          if frmin < bestfl(nbest) & bestla~=1  then
            bestxl(:,nbest) = xrmin;
            bestfl(nbest)   = frmin;
            [bestls,bestis] = sort(-bestfl); bestls = -bestls;
            xltmp  = bestxl(:,bestis); bestxl=xltmp;
            bestfl = bestls;
          end
          for j=1:r; k=j+p; x(:,k)=xr(:,j); fv(k)=fr(j); end
          accept = 1;
        elseif prob >= rand() then
          for j=1:r; k=j+p; x(:,k)=xr(j); fv(k)=fr(j); end
          accept = 2;
        elseif r==n then
          for j=1:r; k=j+p; x(:,k)=xr(j); fv(k)=fr(j); end
          accept = 3;
        else
          r      = r+1;
          accept = 0;
        end
      end
      xtmp    = zeros(n,n+1);
      [fs,is] = sort(-fv); fs = -fs; xtmp=x(:,is); x=xtmp; fv=fs;
      maxdf = fv(n+1)-fv(1);
      rditc = rditc+1;
      itc   = itc+1;         

      if x_history_defined then
        x_history($+1) = list();
        for ii=1:size(x,2)
          x_history($)($+1) = x(:,ii);
        end
      end
    end
    T = T*rd;      
  end
  xnew(:,1) = x1;
  for j=2:n+1; xnew(:,j)=xnew(:,1)+edge(j-1)*id(:,j-1); end
  x = xnew;
  for j=1:n+1; fv(j) = f(x(:,j)); end;
  fcount  = fcount+n+1;
  xtmp    = zeros(n,n+1);
  [fs,is] = sort(-fv); fs = -fs; xtmp=x(:,is); x=xtmp; fv=fs;
  T     = Tmax;
  maxdf = fv(n+1)-fv(1);
  itc   = 0;
  
  if x_history_defined then
    x_history($+1) = list();
    for ii=1:size(xnew,2)
      x_history($)($+1) = xnew(:,ii);
    end
  end
end

finalbestxl(:,1) = bestxl(:,1);
finalbestfl(1)   = bestfl(1);
finalnbest       = 1;

for j=2:nbest
  accept = 0;
  for i=1:j-1;
    if abs(bestxl(:,j)-bestxl(:,i))<0.1 then
      accept = 1;
    end
  end
  if accept==0 then
    finalnbest                = finalnbest+1;
    finalbestxl(:,finalnbest) = bestxl(:,j);
    finalbestfl(finalnbest)   = bestfl(j);
  end   
end

xmin = finalbestxl(:,1);
fmin = bestfl(1);
edge = edgebest*ones(n,1);

for k=1:nbest;
  x1      = zeros(n,n+1);
  x1(:,1) = bestxl(:,k);
  id      = eye(n,n);
  for j=2:n+1; x1(:,j)=x1(:,1)+edge(j-1)*id(:,j-1); end
  [x,fcoun,f1] = optim_dssa_nelder(x1,f,tol);
  if f1 < fmin then
    xmin = x(:,1);
    fmin = f1;
  end
  fcount = fcount+fcoun;
  
  if (x_history_defined) then
    x_history($+1) = list();
    for ii=1:size(x,2)
      x_history($)($+1) = x(:,ii);
    end
  end
end      
endfunction

///////////////////////
// optim_dssa_nelder //
///////////////////////

function [x,fcount,f1,lhist,histout,simpdata] = optim_dssa_nelder(x0,f,tol,maxit,budget)
//
// Nelder-Mead optimizer, No tie-breaking rule other than MATLAB's sort
//
// C. T. Kelley, December 12, 1996
//
//
// This code comes with no guarantee or warranty of any kind.
//
// function [x,lhist,histout,simpdata] = nelder(x0,f,tol,maxit,budget)
//
// inputs:
//	vertices of initial simplex = x0 (n x n+1 matrix)
//          The code will order the vertices for you and no benefit is
//          accrued if you do it yourself.
//
//       objective function = f
//
//       termination tolerance = tol
//       maximum number of iterations = maxit (default = 100)
//           As of today, dist = | best value - worst value | < tol
//           or when maxit iterations have been taken
//       budget = max f evals (default=50*number of variables)
//                 The iteration will terminate after the iteration that
//                 exhausts the budget
//
//
// outputs:
//	final simplex = x (n x n+1) matrix
//
//       number of iterations before termination = itout (optional)
//       iteration histor = histout itout x 5
//         histout = iteration history, updated after each nonlinear iteration
//                 = lhist x 5 array, the rows are
//                   [fcount, fval, norm(grad), dist, diam]
//                   fcount = cumulative function evals
//                   fval = current best function value
//                   norm(grad) = current simplex grad norm
//                   dist = difference between worst and best values
//                   diam = max oriented length
//       simpdata = data for simplex gradient restart 
//              = [norm(grad), cond(v), bar f]

[nargout, nargin] = argn();

//
// initialize counters
//
lhist=0; fcount=0;
//
// set debug=1 to print out iteration stats
//
_debug=0;
//
// Set the N-M parameters
//

if ~isdef('budget','local') then
  budget = [];
end
if ~isdef('maxit','local') then
  maxit = [];
end

rho=1; chi=2; _gamma=.5; sigma=.5;
dsize = size(x0); n=dsize(1);
if nargin < 4 maxit=100*n; end
if nargin < 5 budget=200*n; end
if n >= 10; maxit = 10*maxit; budget=10*budget; end
//
// set the paramters for stagnation detection/fixup
// setting oshrink=0 gives vanilla Nelder-Mead
//
oshrink=1; restartmax=3; restarts=0;
//
//
// Order the vertices for the first time
//
x=x0; [n,m]=size(x); histout=zeros(maxit*3,5); simpdata=zeros(maxit,3);
itout=0; _orth=0;
xtmp=zeros(n,n+1); z=zeros(n,n); delf=zeros(n,1);
for j=1:n+1; fv(j)=f(x(:,j)); end; fcount=fcount+n+1;
[fs,is]=sort(-fv); fs = -fs; xtmp=x(:,is); x=xtmp; fv=fs;
f1=fv(1);
itc=0; dist=fv(n+1)-fv(1);
diam=zeros(n,1);
for j=2:n+1
   v(:,j-1)=-x(:,1)+x(:,j);
   delf(j-1)=fv(j)-fv(1);
   diam(j-1)=norm(v(:,j-1));
end
sgrad=v'\delf; alpha=1.d-4*max(diam)/norm(sgrad);
lhist=lhist+1;
histout(lhist,:)=[fcount, fv(1), norm(sgrad,%inf), 0, max(diam)];
//
// main N-M loop
//
while(itc < maxit & dist > tol & restarts < restartmax & fcount <= budget)
    fbc=sum(fv)/(n+1);
    xbc=sum(x')'/(n+1);
    sgrad=v'\delf;
    simpdata(itc+1,1)=norm(sgrad);
    simpdata(itc+1,2)=cond(v);
    simpdata(itc+1,3)=fbc;
    if(det(v) == 0) then
        disp('simplex collapse')
        break
    end
    happy=0; itc=itc+1; itout=itc;
//
// reflect
//
    y=x(:,1:n);
    xbart = sum(y')/n;  // centriod of better vertices
    xbar=xbart';
    xr=(1 + rho)*xbar - rho*x(:,n+1);
    fr=f(xr); fcount=fcount+1;
    if(fr >= fv(1) & fr < fv(n)) then happy = 1; xn=xr; fn=fr; end;
//    if(happy==1) disp(' reflect '); end
//
// expand
//
    if(happy == 0 & fr < fv(1)) then
        xe = (1 + rho*chi)*xbar - rho*chi*x(:,n+1);
        fe=f(xe); fcount=fcount+1;
        if(fe < fr) then xn=xe;  fn=fe; happy=1; end
        if(fe >=fr) then xn=xr;  fn=fr; happy=1; end
//        if(happy==1) disp(' expand '); end
    end
//
// contract
//
   if(happy == 0 & fr >= fv(n) & fr < fv(n+1)) then
//
// outside contraction
//
       xc=(1 + rho*_gamma)*xbar - rho*_gamma*x(:,n+1);
       fc=f(xc); fcount=fcount+1;
       if(fc <= fr) then xn=xc; fn=fc; happy=1; end;
//       if(happy==1) disp(' outside '); end;
   end
//
// inside contraction
//
   if(happy == 0 & fr >= fv(n+1)) then
       xc=(1 - _gamma)*xbar+_gamma*x(:,n+1);
       fc=f(xc); fcount=fcount+1;
       if(fc < fv(n+1)) then happy=1; xn=xc; fn=fc; end;
//   if(happy==1) disp(' inside '); end;
   end
//
//  test for sufficient decrease, 
//  do an oriented shrink if necessary
//
   if(happy==1 & oshrink==1)
       xt=x; xt(:,n+1)=xn; ft=fv; ft(n+1)=fn;
//       xt=x; xt(:,n+1)=xn; ft=fv; ft(n+1)=f(xn); fcount=fcount+1;
       fbt=sum(ft)/(n+1); delfb=fbt-fbc; armtst=alpha*norm(sgrad)^2;
       if(delfb > -armtst/n) then
           restarts=restarts+1;
           _orth=1; diams=min(diam);
           sx=.5+sign(sgrad); sx=sign(sx);
           if _debug==1 then
               [itc, delfb, armtst]
           end
           happy=0;
           for j=2:n+1
             x(:,j)   = x(:,1); 
             x(j-1,j) = x(j-1,j)-diams*sx(j-1);
           end
       end
   end
//
//  if you have accepted a new point, nuke the old point and
//  resort
//
   if(happy==1) then
       x(:,n+1)=xn; fv(n+1)=fn;
//       x(:,n+1)=xn; fv(n+1)=f(xn); fcount=fcount+1;
       [fs,is]=sort(-fv); fs = -fs; xtmp=x(:,is); x=xtmp; fv=fs;
   end
//
// You're in trouble now! Shrink or restart.
//
   if(restarts >= restartmax) then disp(' stagnation in Nelder-Mead'); end;
   if(happy == 0 & restarts < restartmax) then
       if(_orth ~=1) then disp(' shrink '); end;
       if(_orth ==1) then
       if _debug == 1 then disp(' restart '); end
       _orth=0; end;
       for j=2:n+1;
           x(:,j)=x(:,1)+sigma*(x(:,j)-x(:,1));
           fv(j)=f(x(:,j));
       end
       fcount=fcount+n;
       [fs,is]=sort(-fv); fs = -fs; xtmp=x(:,is); x=xtmp; fv=fs;
   end
//
//  compute the diameter of the new simplex and the iteration data
//
   for j=2:n+1
       v(:,j-1)=-x(:,1)+x(:,j);
       delf(j-1)=fv(j)-fv(1);
       diam(j-1)=norm(v(:,j-1));
   end
   dist=fv(n+1)-fv(1);
   lhist=lhist+1;
   sgrad=v'\delf;
   histout(lhist,:)=[fcount, fv(1), norm(sgrad,%inf), dist, max(diam)];
   f1=fv(1);
end 
endfunction

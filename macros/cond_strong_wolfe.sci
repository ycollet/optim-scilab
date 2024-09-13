function [answer] = cond_strong_wolfe(xk, pk, strong_wolfe_f, strong_wolfe_fs, fk, fsk, alpha)
// Function to test the strong wolfe condition
// xk    : point where to test the strong wolf condition
// pk    : search direction
// strong_wolfe_f     : objective function
// strong_wolfe_fs    : derivative function of the objective function
// fk    : value of the objective function  at point xk
// alpha : value of the step to test

c1 = 1E-4; // c1 in [0,1]
c2 = 0.9;  // 1<c1<c2<1

answer = %T;
// Sufficient decrease condition
answer = answer & (strong_wolfe_f(xk+alpha*pk) <= fk + c1*alpha*fsk'*pk);
// Curvature condition
answer = answer & (abs(strong_wolfe_fs(xk+alpha*pk)'*pk)<=c2*abs(fsk'*pk));
endfunction

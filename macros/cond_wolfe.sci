function [answer] = cond_wolfe(xk, pk, wolfe_f, wolfe_fs, fk, fsk, alpha)
// Function to test the wolfe condition
// xk    : point where to test the strong wolf condition
// pk    : search direction
// wolfe_f     : objective function
// wolfe_fs    : derivative function of the objective function
// fk    : value of the objective function  at point xk
// alpha : value of the step to test

c1 = 1E-4; // c1 in [0,1]
c2 = 0.9;  // 1<c1<c2<1

answer = %T;
// Sufficient decrease condition
answer = answer & (wolfe_f(xk+alpha*pk) <= fk + c1*alpha*fsk'*pk);
// Curvature condition
answer = answer & (wolfe_fs(xk+alpha*pk)'*pk>=c2*fsk'*pk);
endfunction

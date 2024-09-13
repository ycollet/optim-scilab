function [answer] = cond_goldstein(xk, pk, goldstein_f, goldstein_fs, fk, fsk, alpha)
// Function to test the goldstein condition
// xk    : point where to test the wolf condition
// pk    : search direction
// goldstein_f     : objective function
// goldstein_fs    : derivative function of the objective function
// fk    : value of the objective function  at point xk
// alpha : value of the step to test

c = 1E-2; // 0<c<1/2

answer = %T;
f_test = goldstein_f(xk+alpha*pk);

// Sufficient decrease condition
answer = answer & (fk + (1-c)*alpha*fsk'*pk<=f_test);
// Curvature condition
answer = answer & (f_test<=fk+c*alpha*fsk'*pk);
endfunction

function [answer] = cond_armijo(xk, pk, armijo_f, armijo_fs, fk, fsk, alpha)
// Function to test the armijo condition
// xk    : point where to test the armijo condition
// pk    : search direction
// armijo_f     : objective function
// armijo_fs    : derivative function of the objective function
// fk    : value of the objective function  at point xk
// fsk   : value of the derivative function at point xk
// alpha : value of the step to test

c1 = 1E-4; // c1 in [0,1]

answer = %T;
// Sufficient decrease condition (aka Armijo condition)
answer = (armijo_f(xk+alpha*pk) <= fk + c1*alpha*fsk'*pk);
endfunction

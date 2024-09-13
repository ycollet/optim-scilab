function [answer] = cond_gradient(xk, pk, gradient_f, gradient_fs, fk, fsk, alpha)
// Function to test the armijo condition
// xk    : point where to test the armijo condition
// pk    : search direction
// gradient_f     : objective function
// gradient_fs    : derivative function of the objective function
// fk    : value of the objective function  at point xk
// fsk   : value of the derivative function at point xk
// alpha : value of the step to test

c1 = 1E-4;

answer = %T;
// Sufficient null gradient condition
answer = norm(gradient_fs(xk+alpha*pk))<c1;
endfunction

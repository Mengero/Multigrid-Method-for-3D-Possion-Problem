function [A] = tridiag(n)
    e=ones(n,1); h=1/(n+1); h2i=1./(h*h);
    A=h2i*spdiags([-e 2*e -e],-1:1,n,n);
end
function [X] = spline5(u, bound)
    u0 = u(1);
    u1 = u(2);
    b1 = bound(1);
    b2 = bound(2);
    b3 = bound(3);
    b4 = bound(4);
    b5 = bound(5);
    b6 = bound(6);
    
    V = [b1 b2 b3 b4 b5 b6]';
    M = [u0^5 u0^4 u0^3 u0^2 u0 1;
        u1^5 u1^4 u1^3 u1^2 u1 1;
        5*u0^4 4*u0^3 3*u0^2 2*u0 1 0;
        5*u1^4 4*u1^3 3*u1^2 2*u1 1 0;
        20*u0^3 12*u0^2 6*u0 2 0 0;
        20*u1^3 12*u1^2 6*u1 2 0 0];
    X = M\V; % == inv(M)*V
end
    
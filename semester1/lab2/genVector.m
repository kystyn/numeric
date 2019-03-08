function [v, ev, evvec] = genVector( N, ideal )   
    ev = zeros(1, N);
    v = zeros(1, N);
    for i = 1 : 1 : N
        ev(i) = 0.01 * rand() * ~ideal;
        v(i) = 1 + ev(i);
    end
    evvec = ev;
    ev = norm(ev);
end
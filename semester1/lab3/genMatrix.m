function A = genMatrix( N )
    v = zeros(1, N);
    for i = 1 : 1 : N
        v(i) = 50 * rand();
    end
        
    A = full(sprandsym(N, 0.92, v));
end
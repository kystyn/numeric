function v = genVector( N )   
    v = zeros(1, N);
    for i = 1 : 1 : N
        v(i) = rand() * N;
    end
end
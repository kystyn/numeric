function N = prepare(a, b)
F = fopen('lsm.in', 'w');
N = 0;

grids(1) = 'R';
grids(2) = 'U';
grids(3) = 'C';

for gridT = 1 : 1 : 3
    for funcT = 0 : 1 : 1
        for weightT = 0 : 1 : 2
            for nodeCount = 2 : 1 : ~funcT * 20 + funcT * 40
              fprintf(F, '%f %f %i %c %i %i\n', a, b, ...
                nodeCount, grids(gridT), funcT, weightT);
            N = N + 1;
            end
        end
    end
end

fclose(F);
end
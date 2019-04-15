function prepare(F, a, b, funcT, tollerance )
    fprintf(F, '%1.16f %1.16f %i %f\n', a, b, funcT, log10(tollerance));
end
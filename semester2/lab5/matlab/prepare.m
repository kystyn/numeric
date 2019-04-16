function prepare(F, a, b, cp, tollerance )
    fprintf(F, '%1.16f %1.16f %1.16f %i %f\n', a, b, cp, 10, log10(tollerance));
end
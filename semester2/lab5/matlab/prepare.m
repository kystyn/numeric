function prepare(F, a, b, cp, tollerance )
    fprintf(F, '%1.16f %1.16f %1.16f %f\n', a, b, cp, log10(tollerance));
end
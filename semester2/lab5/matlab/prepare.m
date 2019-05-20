function prepare(F, t, a, b, cp, tollerance )
    fprintf(F, '2 %i %1.16f %1.16f %1.16f %1.16f %i %1.16f\n', t, a, b, cp, 10, tollerance);
end
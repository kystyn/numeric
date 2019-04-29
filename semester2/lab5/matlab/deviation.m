function res = deviation( f1, X, f2, N ) %f1 is grid function
res = 0;
for i = 1 : 1 : N
    if (abs(f1(i) - f2(X(i))) > res)
        res = abs(f1(i) - f2(X(i)));
    end
end
end
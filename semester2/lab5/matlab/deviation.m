function res = deviation( f1, f2, a, b, N ) %f1 is grid functiob
res = 0;
h = (b - a) / (N - 1);
for i = 1 : 1 : N
    if (abs(f1(i) - f2(a + (i - 1) * h)) > res)
        res = abs(f1(i) - f2(a + (i - 1) * h));
    end
end
end
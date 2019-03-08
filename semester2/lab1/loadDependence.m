function loadDependence(fname, color)
F = fopen(fname);

format long;
x = fscanf(F, '%f');

n = size(x);

y = x(2 : 2 : n);
x = x(1 : 2 :  n - 1);

semilogy(x, y, color);
title('«ависимость точности от числа узлов');
xlabel('Node count');
ylabel('Deviation');
end
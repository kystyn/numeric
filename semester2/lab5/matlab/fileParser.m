function [X, Y, minfrag, maxfrag, N] = fileParser( fileName )
F = fopen(fileName, 'r');

data = fscanf(F, '%f');

minfrag = data(1);
maxfrag = data(2);

N = length(data);

X = data(3 : 2 : N - 1);
Y = data(4 : 2 : N);

N = (N - 2) / 2;

end
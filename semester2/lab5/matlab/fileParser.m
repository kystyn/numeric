function [X, Y, frag, minfrag, maxfrag, N] = fileParser( fileName )
F = fopen(fileName, 'r');

data = fscanf(F, '%f');

frag = data(1);
minfrag = data(2);
maxfrag = data(3);

N = length(data);

X = data(4 : 2 : N - 1);
Y = data(5 : 2 : N);

N = (N - 2) / 2;

end
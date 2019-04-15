function [X, Y, minfrag, maxfrag] = fileParser( fileName )
F = fopen(fileName, 'r');

data = fscanf(F, '%f');

minfrag = data(1);
maxfrag = data(2);

X = data(3 : 2 : 3 + (maxfrag - 1) * 2);
Y = data(4 : 2 : 4 + (maxfrag - 1) * 2);

end
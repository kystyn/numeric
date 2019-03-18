function [X, Y] = fileParser( fileName )
F = fopen(fileName, 'r');

data = fscanf(F, '%f');

Y = data(1 : 2 : data.length - 1);
X = data(2 : 2 : data.length);

end
function [X, Y] = fileParser( fileName )
F = fopen(fileName, 'r');

data = fscanf(F, '%f');

X = data(1 : 2 : length(data) - 1);
Y = data(2 : 2 : length(data));

end
function [X, Y, frag, minfrag, maxfrag, errlocal, volume, N] = fileParserBoundary( fileName )
F = fopen(fileName, 'r');

data = fscanf(F, '%f');

frag = data(1);
minfrag = data(2);
maxfrag = data(3);
errlocal = data(4);
volume = data(5);

N = length(data);

X = data(6 : 2 : N - 1);
Y = data(7 : 2 : N);

N = (N - 5) / 2;

fclose(F);

end
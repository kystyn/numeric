function [mask, nodeCounts, deviations] = fileParser( fileName, N )
F = fopen(fileName, 'r');

% x = zeros(1, N);
% y = zeros(1, N);
nodeCounts = zeros(1, N);
deviations = zeros(1, N);
mask = zeros(N, 3);

i = 1;
format long;

while (1)
    data = fscanf(F, '%c%i%i', [1 3]);
    if (feof(F))
        break;
    end
    mask(i, 1 : 3) = data;
    steps = fscanf(F, '%i', 1);
    %dots = fscanf(F, '%f', steps * 2);
    %x(i) = dots(1 : steps, 1);
    %y(i) = dots(1 : steps, 2);
    nodeCounts(i) = fscanf(F, '%i', 1);
    deviations(i) = fscanf(F, '%f', 1);
    fscanf(F, '%c', 1);
    i = i + 1;
end

end
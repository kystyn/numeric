function [x, y, deviation] = loadPolinom( fileData )
format long;
n = size(fileData) - 1;

deviation = x(n + 1);

y = x(2 : 2 : n);
x = x(1 : 2 :  n - 1);
end

function [x0, omega, epsilon] = prepare( testCnt )
N = 15;
F = fopen('m.in', 'w');   

v0 = genVector(N, 1);
omega = zeros(1, testCnt);
x0 = zeros(testCnt, N);

% x0
omg = 1.3;
eps = 1e-15;
M = genMatrix(N);
for i = 1 : 1 : testCnt
  omega(i) = 0.1 * i / testCnt + 1.9 * (1 - i/ testCnt);
  fprintf(F, '%i\n', N);
  fprintf(F, '%.16f %.16f\n', omg, eps);
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n');
  fprintf(F, '%5.16f ', transpose(v0));
  fprintf(F, '\n');
  x0(i, 1: N) = genVector(N, i);
  fprintf(F, '%5.16f ', x0(i, 1: N));
  fprintf(F, '\n\n');
end
   
% omega
eps = 1e-10;
for i = 1 : 1 : testCnt
  omega(i) = 0.1 * i / testCnt + 1.9 * (1 - i/ testCnt);
  fprintf(F, '%i\n', N);
  fprintf(F, '%.16f %.16f\n', omega(i), eps);
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n');
  fprintf(F, '%5.16f ', transpose(v0));
  fprintf(F, '\n');
  fprintf(F, '%5.16f ', v0);
  fprintf(F, '\n\n');
end

% epsilon
epsilon = zeros(1, testCnt);
for i = 1 : 1 : testCnt
  epsilon(i) = 10 ^ (-15 * i / testCnt + -2 * (1 - i/ testCnt));
  fprintf(F, '%i\n', N);
  fprintf(F, '%.16f %.16f\n', omg, epsilon(i));
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n');
  fprintf(F, '%5.16f ', transpose(v0));
  fprintf(F, '\n');
  fprintf(F, '%5.16f ', v0);
  fprintf(F, '\n\n');
end
  fprintf(F, '\n');
  fclose(F);

end

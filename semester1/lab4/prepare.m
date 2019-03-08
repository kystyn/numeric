function [v2, epsilon, shifts] = prepare( testCnt )
N = 15;
F = fopen('m.in', 'w');   

shifts = zeros(1, testCnt);
v2 = zeros(1, testCnt);

v = zeros(1, N);
%eigenvalues
for i = 1 : 1 : N
      v(i) = 1 * i;
end
v0 = v;

% separability
eps = 1e-10;
shift = 1.1;
for i = 1 : 1 : testCnt
  v(2) = v0(1) + v0(1) * 10 ^ (-14.0 * (i - 1) / (testCnt - 1));
  v2(i) = v(2);
  M = genMatrix(N, v);
  disp(det(M));
  fprintf(F, '%i\n', N);
  fprintf(F, '%.16f\n', eps);
  fprintf(F, '%.16f\n', shift);
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n\n');
end
 
% epsilon
shift = 0;
epsilon = zeros(1, testCnt);
M = genMatrix(N, v0);
for i = 1 : 1 : testCnt
  epsilon(i) = 10 ^ (-15 * i / testCnt + -2 * (1 - i/ testCnt));
  fprintf(F, '%i\n', N);
  fprintf(F, '%.16f \n', epsilon(i));
  fprintf(F, '%.16f \n', shift);
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n\n');
end

% shift
eps = 1e-10;
for i = 1 : 1 : testCnt
              %
  shifts(i) = (v0(N) * 0.9) * (1 - double(i - 1) / (testCnt - 1)) + (0) * double(i - 1) / (testCnt - 1);
  fprintf(F, '%i\n', N);
  fprintf(F, '%.16f \n', eps);
  fprintf(F, '%.16f \n', shifts(i));
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n\n');
end

fprintf(F, '\n');
fclose(F);

end

function dets = prepare2( testCnt )
N = 15;
F = fopen('m.in', 'w'); 

format long;

dets = zeros(1, N);

omg = 1.2;
eps = 1e-10;
M0 = genMatrix(N);
v0 = genVector(N, 1);
detM0 = det(M0);
for i = 1 : 1 : testCnt
  d = 10 ^ (-0.5 * i);
  fprintf(F, '%i\n', N);
  fprintf(F, '%.16f %.16f\n', omg, eps);
  %M = M0 * ((d / detM0) ^ (1.0 / N));
  M0 = genMatrix(N);
  M = M0 * ((d / det(M0)) ^ (1.0 / N));
  dets(i) = det(M);
  disp(det(M));
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n');
  fprintf(F, '%5.16f ', M * transpose(v0));
  fprintf(F, '\n');
  x0(i, 1: N) = genVector(N, i);
  fprintf(F, '%5.16f ', x0(i, 1 : N));
  fprintf(F, '\n\n');
end
end
function [conds, ev] = prepareSLAEPlotCond()
N = 15;
F = fopen('m.in', 'w');   

v0 = genVector(N, 1);
[v, ev, evvec] = genVector(N, 0);
ev = ev / norm(v0);
disp(ev);
   
for i = 1 : 1 : 1000
  fprintf(F, '%i\n', N);
  M = genRandMatrix(N, 10 ^ (0.013 * i));
  conds(i) = cond(M);
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n');
  fprintf(F, '%5.16f ',  M * transpose(v0));
  %fprintf(F, '%5.16f ',  transpose(v0));
  fprintf(F, '\n\n');
  fprintf(F, '%i\n', N);
  %M = genRandMatrix(N, 10 ^ (0.01 * i));
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n');
  %[a b evvec] = genVector(N, 0);
  fprintf(F, '%5.16f ',  M * (transpose(v0) + transpose(evvec)));
  %fprintf(F, '%5.16f ',  transpose(v));
  fprintf(F, '\n\n');
end
  fprintf(F, '\n');
  fclose(F);
end

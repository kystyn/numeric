function [conds, ev] = prepareSLAEPlotXBerr()
N = 15;
F = fopen('m.in', 'w');   
  
M = genRandMatrix(N, 1e6);
conds = cond(M);

fprintf(F, '%i\n', N);
v0 = genVector(N, 1);
fprintf(F, '%5.16f ', M);
fprintf(F, '\n');
fprintf(F, '%5.16f ', M * transpose(v0));
fprintf(F, '\n');
  
for i = 1 : 1 : 1000
  fprintf(F, '%i\n', N);
  [v, ev(i)] = genVector(N, ~i);
  fprintf(F, '%5.16f ', M);
  fprintf(F, '\n');
  fprintf(F, '%5.16f ', v);
  fprintf(F, '\n');
end
fclose(F);
end
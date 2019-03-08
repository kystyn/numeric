function A = genMatrix( N, v )       
   A = full(sprandsym(N, 0.92, v));
   disp(eig(A));
   %A = D;
end
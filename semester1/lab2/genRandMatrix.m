function M = genRandMatrix( N, cnd )
   M = zeros(N, N);
   A = zeros(N, N);
   
   for i = 1 : 1 : N
      M(i, i) = 1;
      A(i, i) = 1;
   end
   M(1, 1)  = cnd;
    
   r = rand(N, 1, 'double');
   r = r / norm(r);
   A = A - 2 * r * transpose(r);
   M = A * M * transpose(A) / cnd;
end

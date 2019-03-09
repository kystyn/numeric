function main()
 N = prepare(0.1, 1.1);
 system('./lab1');
 [mask, nodeCounts, deviations] = fileParser('hermit.out', N);
 
 chkMask = zeros(1, 3);
 color(1) = 'r';
 color(2) = 'b';
 color(3) = 'g';
 
figure;
 
grids(1) = 'R';
grids(2) = 'U';
grids(3) = 'C';

 for i = 1 : 1 : 3
     for j = 0 : 1 : 1
         chkMask(1) = grids(i);
         chkMask(2) = 1;
         chkMask(3) = j;
         clr = color(i);
         if (j == 0)
             clr = strcat(clr, '-');
         end
         drawDependencies(mask, nodeCounts, deviations, chkMask, clr);
         hold on;
     end
 end
 
 legend('Random, non-diff', 'Random, diff', ...
        'Uniform, non-diff', 'Uniform, diff', ...
        'Chebyshev, non-diff', 'Chebyshev, diff');

end
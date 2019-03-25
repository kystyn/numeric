function run()
  format long;
  a = [0, -pi / 3];
  b = [pi / 2, pi / 3];
  
  fact = 1;
  
  N = 8;
  tolls = zeros(1, N);
  facts = tolls;
  for funcT = 0 : 1 : 1
    F = fopen('integral.in', 'w');
    for toll = 1 : 1 : N
      prepare(F, a(funcT + 1), b(funcT + 1), funcT, 10 ^ -toll);
      tolls(toll) = 10 ^ -toll;
      facts(toll) = fact;
    end
    fclose(F);
    
    %system('lab3.exe');
    [X, Y] = fileParser('integral.out');
    
    figure;
    drawSmth(tolls, X - facts, 'r', strcat('Dependence of fact precision on settable ', char('0' + funcT), '-diff. function'), ...
    'Settable precision', 'Fact precision');
    grid on;

    figure;
    drawSmth(X - facts, Y, 'r', strcat('Dependence of stpe count on precision', char('0' + funcT),'-diff. function'), ...
    'Fact precision', 'Step count');
    grid on;
  end
end
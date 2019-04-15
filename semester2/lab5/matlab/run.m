function run()
  format long;
  a = [0, -pi / 3];
  b = [pi / 2, pi / 3];
    
  N = 8;
  tolls = zeros(1, N);
  fact = zeros(1, N);
  for funcT = 0 : 1 : 1
    F = fopen('integral.in', 'w');
    for toll = 1 : 1 : N
      prepare(F, a(funcT + 1), b(funcT + 1), 1 - funcT, 10 ^ -toll);
      tolls(toll) = 10 ^ (0 - toll);
      fact(toll)=1;
    end
    fclose(F);
    
    %system('lab3.exe');
    system('./lab3');
    [X, Y] = fileParser('integral.out');
    
    figure;
    drawSmth(tolls, abs(X - fact), 'r', strcat('Dependence of fact precision on settable ', char('0' + funcT), '-diff. function'), ...
    'Settable precision', 'Fact precision');
    grid on;

    figure;
    drawSmth(tolls, Y, 'r', strcat('Dependence of step count on precision ', char('0' + funcT),'-diff. function'), ...
    'Fact precision', 'Step count');
    grid on;
  end
end
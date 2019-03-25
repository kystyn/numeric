function run()
  format long;
  a = [0, -pi / 3];
  b = [pi / 2, pi / 3];
  
  fact = 1;
  
  tolls = zeros(1, 16);
    for funcT = 0 : 1 : 1
    F = fopen('integral.in', 'w');
    for toll = 1 : 1 : 16
      prepare(F, a(funcT + 1), b(funcT + 1), funcT, 10 ^ -toll);
      tolls(toll) = 10 ^ -toll;
    end
    fclose(F);
    
    system('./lab3');
    [X, Y] = fileParser('integral.out');
    
    grid on;
    
    figure;
    drawSmth(tolls, X - tolls, 'r', strcat('Dependence of fact precision on settable ', char('0' + funcT),'-diff. function'));
    figure;
    drawSmth(X - tolls, Y, strcat('Dependence of stpe count on precision', char('0' + funcT),'-diff. function'));
  end
end
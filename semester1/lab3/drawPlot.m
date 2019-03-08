function drawPlot( testCnt, x0, omega, epsilon )
    F = fopen('solve.out', 'r');         
    x = fscanf(F, '%f');
    fclose(F);
    
    X = zeros(1, testCnt);
    Y = zeros(1, testCnt);
    
    %----------
    cnt = 1;
    for li = 1 : 1 : testCnt
        N = x(cnt);
        root = x(cnt + 1 : cnt + N);
        cnt = cnt + N + 1;
        st = x(cnt);
        cnt = cnt + 1;
        
        X(li) = norm(root - transpose(x0(li, 1 : N)));
        Y(li) = st;
    end
    figure;
    plot(x0, Y, 'r');
    title('Dependence of steps count on x_0');
    xlabel('||x_0-x*||');
    ylabel('Steps');
    
    %----------
    for li = 1 : 1 : testCnt
        N = x(cnt);
        cnt = cnt + N + 1;
        st = x(cnt);
        cnt = cnt + 1;

        X(li) = omega(li);
        Y(li) = st;
    end
    figure;
    plot(omega, Y, 'r');
    title('Dependence of steps count on \omega');
    xlabel('\omega');
    ylabel('Steps');
    
    %----------
    for li = 1 : 1 : testCnt
        N = x(cnt);
        cnt = cnt + N + 1;
        st = x(cnt);
        cnt = cnt + 1;

        Y(li) = st;
    end
    figure;
    semilogx(epsilon, Y, 'r');
    title('Dependence of steps count on \epsilon');
    xlabel('\epsilon');
    ylabel('Steps');
end  

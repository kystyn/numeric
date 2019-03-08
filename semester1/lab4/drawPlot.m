function drawPlot( testCnt, v2, epsilon, shift )
    F = fopen('solve.out', 'r');         
    x = fscanf(F, '%f');
    fclose(F);
    
    X = zeros(1, testCnt);
    Y = zeros(1, testCnt);
    accuracy = zeros(1, testCnt);
    
    % separability: steps and true accuracy
    cnt = 1;
    for li = 1 : 1 : testCnt
        X(li) = v2(li) - 1;

        accuracy(li) = abs(x(cnt) - 1);
        cnt = cnt + 1;
        Y(li) = x(cnt);
        cnt = cnt + 1;
    end
    figure;
    semilogx(X, Y, 'r');
    title('Dependence of steps count on separability');
    xlabel('|\lambda_2 - \lambda_1|');
    ylabel('Steps');
    
    figure;
    loglog(X, accuracy, 'r');
    title('Dependence of true accuracy on separability');
    xlabel('|\lambda_2 - \lambda_1|');
    ylabel('|\lambda_1^*-\lambda_1^0|');
            
    % epsilon
    for li = 1 : 1 : testCnt
        cnt = cnt + 1;
        st = x(cnt);
        cnt = cnt + 1;

        Y(li) = st;
    end
    figure;
    semilogx(epsilon, Y, 'r');
    title('Dependence of steps count on \epsilon');
    xlabel('\epsilon');
    ylabel('Steps');
    
    % shift
    for li = 1 : 1 : testCnt
        cnt = cnt + 1;
        st = x(cnt);
        cnt = cnt + 1;

        Y(li) = st;
    end
    figure;
    plot(shift, Y, 'r');
    title('Dependence of steps count on shift');
    xlabel('Shift');
    ylabel('Steps');
end  

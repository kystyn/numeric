function drawPlot2( testCnt, X )
    F = fopen('solve.out', 'r');         
    x = fscanf(F, '%f');
    fclose(F);
    
    Y = zeros(1, testCnt);
    
    %----------
    cnt = 1;
    for li = 1 : 1 : testCnt
        N = x(cnt);
        root = x(cnt + 1 : cnt + N);
        cnt = cnt + N + 1;
        st = x(cnt);
        cnt = cnt + 1;
        
        Y(li) = st;
    end
    figure;
    semilogx(X, Y, 'r');
    title('Dependence of steps count on det(A)');
    xlabel('det(A)');
    ylabel('Steps');
end  

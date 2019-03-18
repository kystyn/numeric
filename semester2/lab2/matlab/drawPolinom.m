function drawPolinom(X, Y, trueY, style, trueStyle)
hold on;
grid off;

plot(X, Y, style);
plot(X, trueY, trueStyle);
legend('Interpolation', 'Original');
title('Interpolated function and original');
xlabel('x');
ylabel('y');
end
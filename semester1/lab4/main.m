function main()
format long;
rng('shuffle');
testCnt = 15;
[v2, epsilon, shift] = prepare(testCnt);
system('lab4.exe');
drawPlot(testCnt, v2, epsilon, shift);
end
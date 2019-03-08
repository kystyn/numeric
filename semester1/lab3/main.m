function main()
format long;
testCnt = 20;
[x0, omega, epsilon] = prepare(testCnt);
system('lab3.exe');
drawPlot(testCnt, x0, omega, epsilon);

[dets] = prepare2(testCnt);
system('lab3.exe');
drawPlot2(testCnt, dets);
end
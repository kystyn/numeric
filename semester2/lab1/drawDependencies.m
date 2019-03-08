function drawDependencies( mask, nodeCount, deviation, chkMask, color )
j = 1;
for i = 1 : 1 : length(nodeCount)
  if (mask(i, 1 : 3) == chkMask(1 : 3))
      X(j) = nodeCount(i);
      Y(j) = deviation(i);
      j = j + 1;
  end
end

drawSmth(X, Y, color, 'Dependence of deviation on node count');

end
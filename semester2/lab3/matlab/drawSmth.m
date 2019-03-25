function drawSmth( X, Y, style, titl )
  semilogy(X, Y, style);
  title(titl);
  xlabel('Node count');
  ylabel('Deviation');
end
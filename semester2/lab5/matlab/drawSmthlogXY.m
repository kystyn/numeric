function drawSmthlogXY( X, Y, style, titl, xl, yl )
  semilogx(X, Y, style);
  title(titl);
  xlabel(xl);
  ylabel(yl);
end
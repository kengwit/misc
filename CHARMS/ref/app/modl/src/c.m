for i=0:5
  for j=i:5
    if (i == j)
       fprintf(1,'C%d%d=C[%d][%d];\n', i,j,i,j);
    else
       fprintf(1,'C%d%d=C%d%d=C[%d][%d];\n', i,j,j,i,i,j) ;
    end
  end
end

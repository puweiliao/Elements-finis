function val = sigma (X,zone)
  x=X(1);
  y=X(2);
  if zone==1
      val=sigma_1(x,y);
    elseif zone==2
      val=sigma_2(x,y);
     else
      val="error";
end

endfunction

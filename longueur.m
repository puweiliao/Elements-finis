function [c] = longueur(xy)         %calcule du longueur
    global A; global B; global L;
    nb = length(L);
    nn = length(xy)/2;      
    x = xy(1:nn);
    y = xy(nn+1:end);
    c(1)=(x(1))^2+(y(1))^2-L(1)^2;
    for i=2:nn
        c(i)=(x(i)-x(i-1))^2+(y(i)-y(i-1))^2-L(i)^2;
    end
    c(nb)=(A-x(nn))^2+(B-y(nn))^2-L(nb)^2;
    c = c';
    return
end

function [df] = DF(phi, i, xy)         %calcule de la diff��rence finie
    n = length(xy);
    t_i = sqrt(eps) * max(1, abs(xy(i)));
    e_i = double(1:n == i);
 
    df = (phi(xy + (t_i*e_i)') - phi(xy - (t_i*e_i)') ) / (2 * t_i);
    return
end




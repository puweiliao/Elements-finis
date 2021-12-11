[e,c,g,a,hl,indic] = chs(4,xy,[0.01,0.01,0.01,0.01,0.01]);
function [err] = erreur_grad(phi, i, xy,g)
    n = length(xy);
    ti = sqrt(eps) * max(1, abs(xy(i)));
    ei = double(1:n == i);
    DF = (phi(xy + (ti*ei)') - phi(xy - (ti*ei)') ) / (2 * ti);
    if g(i) == 0
        err = abs(DF - g(i));
    end
    if g(i) ~= 0
        err = abs(DF - g(i))/g(i);
    end
    printf("\n %d ; pas: %e; f'(x): %e; DF: %e \n;erreur: %e ", i, ti,g(i), DF,err);
end
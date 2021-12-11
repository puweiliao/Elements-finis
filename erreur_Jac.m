function [Jac] = erreur_Jac(phi, xy)
    Jac = sparse(length(phi(xy)), length(xy));
    err_J = sparse(length(phi(xy)), length(xy));
    for i=1:length(xy)
        [Jac(1:end, i)] = DF(phi, i, xy);
    end
    return
end 





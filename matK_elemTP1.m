function [Kel] = matK_elemTP1(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;


% calcul de la matrice de raideur
% -------------------------------
Y=[0,y1-y2,y1-y3; y2-y1,0,y2-y3;y3-y1,y3-y2,0];      

X=[0,x1-x2,x1-x3; x2-x1,0,x2-x3; x3-x1,x3-x2,0];


Kel = zeros(3,3);

Kel(1,2) = 1/(2*D)*(Y(2,3)*Y(3,1)+X(2,3)*X(3,1));
Kel(1,3)=1/(2*D)*(Y(2,3)*Y(1,2)+X(2,3)*X(1,2));
Kel(2,3)=1/(2*D)*(Y(3,1)*Y(1,2)+X(3,1)*X(1,2));

Kel=Kel+Kel.';   %on ajoute la transposée pour éviter de recopier

Kel(1,1)=1/(2*D)*(Y(2,3)^2+X(2,3)^2);
Kel(2,2)=1/(2*D)*(Y(3,1)^2+X(3,1)^2);
Kel(3,3)=1/(2*D)*(Y(1,2)^2+X(1,2)^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020

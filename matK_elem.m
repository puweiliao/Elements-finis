function [Kel] = matK_elem(S1, S2, S3,zone)
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
% NOTE (1) Utilisation d une quadrature a 3 point d ordre 2
%    
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

%quantités qui nous serviront dans la suite

d1=(x2-x1)*(y3-y1);
d2=(y2-y1)*(x3-x1);

B=[[x2-x1,x3-x1];[y2-y1,y3-y1]] ; 

P=inv(transpose(B));

nabla_w=[[-1;-1],[1;0],[0;1]];

S_q=[[1/6;1/3], [1/6;1/3],[1/3;1/6]];

S=S1.';

phi=[B*S_q(:,1)+S , B*S_q(:,2)+S , B*S_q(:,3)+S];


% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i=1:3
    for j=1:3
    for k=1:3
      Kel(i,j) += sigma(phi(:,k),zone);
     end; %k
     Kel(i,j)=abs(D)/6*dot(P*nabla_w(:,i),P*nabla_w(:,j)) * Kel(i,j);     #je ne sais pas pourquoi il faut mettre un -...
  end; % j
end; % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020

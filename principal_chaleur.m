% =====================================================
%
% principal_chaleur;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour 
% 1) l'equation de la chaleur suivante stationnaire, avec condition de
% Dirichlet non homogene
%
% | \alpha T - div(\sigma \grad T)= S,   dans \Omega=\Omega_1 U \Omega_2
% |         T = T_\Gamma,   sur le bord
%
% ou S est la source de chaleur, T_\Gamma la temperature exterieure
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% 2) l'equation de la chaleur dependant du temps avec condition de 
% Dirichlet non homogene
%
% | dT/dt - div(\sigma \grad T)= S,   dans \Omega=\Omega_1 U \Omega_2 et pour tout t< t_max
% |         T = T_\Gamma,   sur le bord et pour tout t< t_max
% |         T = T_0       dans \Omega et pour t=0  
%
% ou S est la source de chaleur, T_\Gamma la temperature exterieure,
% T_0 est la valeur initiale de la temp?rature
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
% =====================================================
% Donnees du probleme
% ---------------------------------
h = 0.05;
system(['gmsh -2 -clmax ' num2str(h) ' -clmin ' num2str(h) ' geomChaleur.geo']);
nom_maillage = 'geomChaleur.msh' ;

validation = 'oui';
pb_stationnaire = 'non';
pb_temporel = 'non';

if strcmp(validation,'oui')
    alpha = 1;
    T_Gamma = 290;
end

if strcmp(pb_stationnaire,'oui')
    alpha = 1;
    T_Gamma = 290;
end

if strcmp(pb_temporel,'oui')
    Tps_initial = 0;
    Tps_final = 1;
    delta_t = 0.01;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    T_Gamma = 290;
end

% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for i=1:Nbtri

  % calcul des matrices elementaires du triangle l 
   k = Numtri(i,1);
   l = Numtri(i,2);
   m = Numtri(i,3);
  
  
   S1=[Coorneu(k,1),Coorneu(k,2)];
   S2=[Coorneu(l,1),Coorneu(l,2)];
   S3=[Coorneu(m,1),Coorneu(m,2)];

   [Kel]=matK_elem(S1,S2,S3,Reftri(i));

   [Mel]=matM_elem(S1,S2,S3);
   
   %On assemble 
   
   MM(k,k) += Mel(1,1);
   MM(l,l) += Mel(2,2);
   MM(m,m) += Mel(3,3);
  
   MM(k,l)+= Mel(1,2);
   MM(k,m)+=Mel(1,3);
  
   MM(l,m)+=Mel(2,3);
   MM(l,k) += Mel(2,1);
  
   MM(m,k) += Mel(3,1);
   MM(m,l)+= Mel(3,2);
  
   %On fait la même chose pour la matrice K.
   KK(k,k) += Kel(1,1);
   KK(l,l) += Kel(2,2);
   KK(m,m) += Kel(3,3);
  
   KK(k,l)+= Kel(1,2);
   KK(k,m)+=Kel(1,3);
  
   KK(l,m)+=Kel(2,3);
   KK(l,k) += Kel(2,1);
  
   KK(m,k) += Kel(3,1);
   KK(m,l)+= Kel(3,2);
   
end % for i

% Matrice EF
% -------------------------
AA = KK;

% =====================================================
% =====================================================
% Pour le probleme stationnaire et la validation
% ---------------------------------

% Calcul du second membre F
% -------------------------
FF=zeros(Nbpt,1);
for i = 1:Nbpt
  S=[Coorneu(i,1),Coorneu(i,2)];
  FF(i) = f(S(1),S(2));
end
LL = MM*FF;


% inversion
% ----------
% tilde_AA ET tilde_LL SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
% APRES PSEUDO_ELIMINATION 

[tilde_AA, tilde_LL] = elimine(AA,LL,Refneu);

T_Gamma=290;
UU = tilde_AA\tilde_LL;
TT = T_Gamma + UU; 

% validation
% ----------
if strcmp(validation,'oui')
   UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
   %affiche(UU_exact, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));

	% Calcul de l erreur L2
	% A COMPLETER
	% Calcul de l erreur H1
	% A COMPLETER
	% attention de bien changer le terme source (dans FF)
end


% visualisation
% -------------
if ( strcmp(validation,'oui') || strcmp(pb_stationnaire,'oui') )
    affiche(TT, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));
end
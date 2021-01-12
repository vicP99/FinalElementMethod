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

h=0.1;
system(['gmsh -2 -clmax ' num2str(h) ' -clmin ' num2str(h) ' geomChaleur.geo']);
nom_maillage = 'geomChaleur.msh';

% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh('geomChaleur.msh');

validation_dirichlet = 'non';
validation_neumann = 'non';
pb_stationnaire_dirichlet = 'non';
pb_stationnaire_neumann = 'non';
pb_temporel = 'non';
pb_temporel_verification = 'non';
pb_temporel_fourier = 'non';
pb_temporel_fourier_verif = 'oui';

if strcmp(validation_dirichlet,'oui')
    alpha = 1;
    T_Gamma = 0;
    lambda=0;  %pour pas considérer le terme des intégrales sur le bord 

    sigma_1=@(vec) 1;
    sigma_2=@(vec) 1;
    f=@(x,y) (1+2*pi^2)*sin(pi*x)*sin(pi*y); 

    U_c=zeros(Nbpt,1); %pour pas considérer le terme des intégrales sur le bord (on fera toujours ca pour dirichlet et neumann)
end

if strcmp(pb_stationnaire_dirichlet,'oui')
    alpha = 1;
    T_Gamma = 290;
    lambda=0;

    sigma_1=@(vec) 5;
    sigma_2=@(vec) sqrt(3)/2;
    f=@(x,y) 600*exp(-((x-1)/0.8).^2-((y-1)/0.8).^2) - T_Gamma; 

    U_c=zeros(Nbpt,1); 
end

if strcmp(validation_neumann,'oui')
    alpha = 1;
    T_Gamma = 0;
    lambda = 1;

    T_c = @(x,y) cos(pi*x)*cos(2*pi*y);
    u_c = @(x,y) T_c(x,y)-T_Gamma;
    U_c = arrayfun(u_c,Coorneu(:,1),Coorneu(:,2));
    
    f=@(x,y) (1+5*pi^2)*cos(pi*x)*cos(2*pi*y);  
    sigma_1=@(vec) 1;
    sigma_2=@(vec) 1;
end

if strcmp(pb_stationnaire_neumann,'oui')
    alpha = 1;
    T_Gamma = 290;
    T_c = 290;
    U_c = (T_c-T_Gamma)*ones(Nbpt,1);
    lambda = 1;

    f=@(x,y) 600*exp(-((x-1)/0.8).^2-((y-1)/0.8).^2) - T_Gamma;  
    sigma_1=@(vec) 5;
    sigma_2=@(vec) 1/4*(2+sin(16*pi*vec(1)))*(2+sin(16*pi*vec(2)));
end

if strcmp(pb_temporel,'oui')
    Tps_initial = 0;
    Tps_final = 1;
    delta_t = 0.01;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    T_Gamma = 280;
    lambda=0;
    T_c=@(x,y,t) 0;

    f=@(x,y,t) 600*exp(-((x-1)/0.8).^2-((y-1)/0.8).^2)*exp(-5*t);  
    sigma_1=@(vec) 5;
    sigma_2=@(vec) 1/4*(2+sin(16*pi*vec(1)))*(2+sin(16*pi*vec(2)));
    condition_initiale=@(x,y) 300;
end
if strcmp(pb_temporel_verification,'oui')
    Tps_initial = 0;
    Tps_final = 11;
    delta_t = 1;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    T_Gamma = 0;
    lambda=0;

    f=@(x,y,t) 5*pi^2*sin(pi*x)*sin(pi*y)*exp(pi^2*t) ;  
    sigma_1=@(vec) 1;
    sigma_2=@(vec) 1;
    condition_initiale=@(x,y) sin(pi*x)*sin(pi*y);
end
if strcmp(pb_temporel_fourier,'oui')
    Tps_initial = 0;
    Tps_final = 1;
    delta_t = 0.01;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    T_Gamma = 280;
    lambda=1;
    T_c = @(x,y,t) 280 -2000*min(t,0.1);

    for i=1:Nbpt
        if Coorneu(i,1) == 0
            Refneu(i)=-1;
        end
    end

    f=@(x,y,t) 600*exp(-((x-1)/0.8).^2-((y-1)/0.8).^2);  
    sigma_1=@(vec) 5;
    sigma_2=@(vec) 1/4*(2+sin(16*pi*vec(1)))*(2+sin(16*pi*vec(2)));
    condition_initiale=@(x,y) 300;
end
if strcmp(pb_temporel_fourier_verif,'oui')
    Tps_initial = 0;
    Tps_final = 1;
    delta_t = 0.1;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    T_Gamma = 0;
    lambda=1;
    T_c = @(x,y,t) cos(9/4*pi*x)*sin(2*pi*y)*sin(4*pi*t*y);
    solReel=@(x,y,t) cos(9/4*pi*x)*sin(2*pi*y)*sin(4*pi*t*y);
    for i=1:Nbpt
        if Coorneu(i,1) == 0
            Refneu(i)=-1;
        end
    end

    f=@(x,y,t) (4+(9/4)^2 + (4*t)^2)*pi^2*cos(9/4*pi*x)*sin(2*pi*y)*sin(4*pi*t*y) + 4*pi*y*cos(9/4*pi*x)*sin(2*pi*y)*cos(pi*t*y) - 16*pi^2*t*cos(9/4*pi*x)*cos(2*pi*y)*cos(4*pi*t*y);  
    sigma_1=@(vec) 1;
    sigma_2=@(vec) 1;
    condition_initiale=@(x,y) 0;
end

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de masse
SS = sparse(Nbpt,Nbpt); % matrice de masse de bord
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  
  % calcul des matrices elementaires du triangle l 
   Mel=matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:));
  
   Kel=matK_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:),sigma_1,sigma_2,Reftri(l));
   
    
    % On fait l'assemblage de la matrice globale
    for i=1:3
        for j=1:3
          MM(Numtri(l,i),Numtri(l,j))= MM(Numtri(l,i),Numtri(l,j))+Mel(i,j);
          KK(Numtri(l,i),Numtri(l,j))=KK(Numtri(l,i),Numtri(l,j))+Kel(i,j);
        end
    end
end % for l

for a=1:Nbaretes
    Sel = mat_elem_surface(Coorneu(Numaretes(a,1),:),Coorneu(Numaretes(a,2),:),Refaretes(a));
    % On fait l'assemblage de la matrice globale
    for i=1:2
        for j=1:2
            SS(Numaretes(a,i),Numaretes(a,j))=SS(Numaretes(a,i),Numaretes(a,j))+Sel(i,j);
        end
    end
end % for a

% Matrice EF
% ------------------------- 
AA = alpha*MM+KK+lambda*SS;



if not(strcmp(pb_temporel,'oui') || strcmp(pb_temporel_verification,'oui')|| strcmp(pb_temporel_fourier,'oui') ||strcmp(pb_temporel_fourier_verif,'oui'))
    % =====================================================
    % =====================================================
    % Pour le probleme stationnaire et la validation
    % ---------------------------------

    % Calcul du second membre F
    % -------------------------
    % A COMPLETER EN UTILISANT LA ROUTINE f.m
    

    FF = arrayfun(f,Coorneu(:,1),Coorneu(:,2));

    LL = MM*FF+lambda*SS*U_c;

    % inversion
    % ----------

    [tilde_AA,tilde_LL] = elimine(AA,LL,Refneu,Nbpt);
    tilde_UU = tilde_AA\tilde_LL;
    tilde_TT = T_Gamma + tilde_UU;

    UU = AA\LL;
    TT =T_Gamma + UU ;
end

% validation
% ----------
if strcmp(validation_dirichlet,'oui')
    UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
	% Calcul de l erreur L2
	(UU_exact-tilde_UU)'*MM*(UU_exact-tilde_UU)
	% Calcul de l erreur H1
	(UU_exact-tilde_UU)'*KK*(UU_exact-tilde_UU)
	% attention de bien changer le terme source (dans FF)
end

if strcmp(validation_neumann,'oui')
    UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
	% Calcul de l erreur L2
	L2=(UU_exact-UU)'*MM*(UU_exact-UU)
	% Calcul de l erreur H1
	H1=(UU_exact-UU)'*KK*(UU_exact-UU)
	% attention de bien changer le terme source (dans FF)
end

% visualisation
% -------------

if strcmp(pb_stationnaire_dirichlet,'oui') || strcmp(validation_dirichlet,'oui')
    affiche(tilde_TT, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));
    hold on;
end


 if ( strcmp(validation_neumann,'oui') || strcmp(pb_stationnaire_neumann,'oui') ) 
    affiche(TT, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
    hold on;
end 



% =====================================================
% =====================================================
% Pour le probleme temporel
% ---------------------------------
if strcmp(pb_temporel,'oui') || strcmp(pb_temporel_verification,'oui') ||strcmp(pb_temporel_fourier,'oui') || strcmp(pb_temporel_fourier_verif,'oui')

    % on initialise la condition initiale
    % -----------------------------------
    T_initial=zeros(Nbpt,1);
    for i=1:Nbpt
        T_initial(i,1) = condition_initiale(Coorneu(i,1),Coorneu(i,2));
    end
	% solution a t=0
    % --------------
    UU = T_initial - T_Gamma ;
    TT = UU+T_Gamma;

    % visualisation
    % -------------
  figure;
    hold on;
    affiche(TT, Numtri, Coorneu, ['Temps = ', num2str(0)]);
    axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2))]);

hold off;

	% Boucle sur les pas de temps
	% ---------------------------
    for k = 1:N_t
        LL_k = zeros(Nbpt,1);
        
        
        % Calcul du second membre F a l instant k*delta t
        % -----------------------------------------------
        % A COMPLETER EN UTILISANT LA ROUTINE f_t.m et le terme precedent (donne par UU)
        F = @(x,y) f(x,y,k*delta_t);
        U_c = @(x,y) T_c(x,y,k*delta_t)-T_Gamma;

        UU_c = arrayfun(U_c,Coorneu(:,1),Coorneu(:,2));
        FF = arrayfun(F,Coorneu(:,1),Coorneu(:,2));


        LL_k = MM*(FF + alpha*UU) + lambda*SS*UU_c;
		% inversion
		% ----------
		% tilde_AA ET tilde_LL_k SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
        % APRES PSEUDO_ELIMINATION 
		% ECRIRE LA ROUTINE elimine.m ET INSERER L APPEL A CETTE ROUTINE
        % A UN ENDROIT APPROPRIE
        
        [tilde_AA,tilde_LL_k]=elimine(AA,LL_k,Refneu,Nbpt);


        UU = tilde_AA\tilde_LL_k;
        TT = UU + T_Gamma;

        if strcmp(pb_temporel_fourier_verif,'oui')
            if abs(k*delta_t - 1) < delta_t 
                Sol = @(x,y) solReel(x,y,k*delta_t);
                TTreel =  arrayfun(Sol,Coorneu(:,1),Coorneu(:,2));
                L2 = (TTreel - TT)'*MM*(TTreel-TT);
                affiche(-TT, Numtri, Coorneu, ['Temps = ', num2str(k*delta_t)]);
                axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2))]); 
                
            end
        end

        % visualisation 
		%& -------------
     pause(0.05)
    affiche(TT, Numtri, Coorneu, ['Temps = ', num2str(k*delta_t)]);
    axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2))]); 
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020


function [Kel] = matK_elem(S1, S2, S3,sigma_1,sigma_2, Reftriangle)
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


  % D est, au signe pres, deux fois l'aire du triangle
  D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
  if (abs(D) <= eps) 
    error('l aire d un triangle est nulle!!!'); 
  end;


  % calcul de la matrice de raideur
  % -------------------------------
  Bl = [x2-x1,x3-x1;y2-y1,y3-y1];
  F_trans=@(x,y) Bl*[x,y]' + [x1,y1]';


  if(Reftriangle==1)
    sigChap=@(x,y) sigma_1(F_trans(x,y));
    else
      sigChap=@(x,y) sigma_2(F_trans(x,y));
    
  end


  function res= grad_w(i)
    if i==1
      res=[-1,-1]';
    end
    if i==2
      res=[1,0]';
    end
    if i==3
      res=[0,1]';
    end
  end


  Kel = zeros(3,3);
  for i=1:3
    for j=1:3
      G=@(x,y) sigChap(x,y)*((Bl'\grad_w(i))'*(Bl'\grad_w(j)))*abs(D);
      Kel(i,j) = 1/6*(G(1/6,1/6)+G(2/3,1/6)+G(1/6,2/3));
    end; % j
  end; % i

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2020

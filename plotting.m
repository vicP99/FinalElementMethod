function [Error,H] = plotting()
    H=0.05:0.01:1;
    Error=zeros(2,size(H,2));
    for i=1:size(H,2)
 
        [L2bis]=principal_chaleur(H(1,i))


        Error(1,i) =(sqrt(L2bis));

%{
         Error(2,i) =(sqrt(H1));  
%}


    end


end
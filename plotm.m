function Error = plot()
    H=0.01:0.01:0.2;
    Error=zeros(2,size(H,2));
    for i=1:size(H,2)
        i
        H(1,i)
        [L2,H1]=principal_chaleur(H(1,i))


        Error(1,i) =log(sqrt(L2));
        Error(2,i) =log(sqrt(H1)); 


    end


end
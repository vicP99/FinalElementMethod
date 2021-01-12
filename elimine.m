function [A,L] = elimine(AA,LL,refneu,Nbpt)

A=AA;
L=LL;
    for i=1:Nbpt
        if refneu(i)==1
            L(i)=0;
            A(i,i)=1;
            for j=1:Nbpt
                if j~=i
                    A(i,j)=0;
                    A(j,i)=0;
                end
            end
        end
    end
end
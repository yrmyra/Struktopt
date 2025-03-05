%% Nummerisk Diff

OPTPAR=optimset('TolX',1e-18,'Display','off');

numdiff_mat=zeros(nelm,1);
for it=1:1
    df=zeros(nelm,1);
%     h=sum(z)/1000000;
    h=1e-6;

    for el=1:nelm
        Kdf=K;
        indx=edof(el,2:end);
        Kdf(indx,indx)=K(indx,indx)+Ki{el}*h;
        adf=solveq(Kdf,F,bc);
        dz=z+z(el)*h;
        df(el)=(g0(adf,F,dz)-g0(a,F,z))/h;
    end
    numdiff_mat(:,it)=df;

    
end
function [dsigma] = sigma_der(Ki,z,u,p,vM,p_sig,P)


dPNdvMz=(sum(vM.^p_sig))^(1/p_sig-1)*vM.^(p_sig-1);
dvMzdvM=sqrt(z);
dvMdsig=1/vM*P*es


dgdu=dPNdvMz*dvMzdvM*dvMdsig*dsigdz;
dKdz=p*(1-delta_0)*z(el)^(p-1)*Ki{el};
lambda=-inv(K)*dgdu;
dsigma=dgdz+lambda'*dKdz*u;

end


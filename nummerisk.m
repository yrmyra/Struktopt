%% 

clc

L=0.01;
t=0.01;
E=4.5*10^9;
nu=0.32; 
V_box=12*L*2*L*t;
meshfac=2^0;
p=3;
last=-800;
maxit=1;
r=0.1*L;



[coord, dof, enod, edof, Ex, Ey, ~] = designDomain(6*L, 2*L, L/meshfac);

nelm=length(edof);
nnod=length(coord);
ndof=nnod*2;

ep=[1 t 2];
D=hooke(1,E,nu);


% hitta noder med bc

last_nod=find(ismember(coord, [0.03, 0.02], 'rows'))*2;


if meshfac>1
    bc_nodes=find(ismember(coord, [0.005, 0;0.055, 0], 'rows'))*2;
else 
    bc_nodes=find(ismember(coord, [0.01, 0;0.05, 0], 'rows'))*2;
end

bc_node=bc_nodes;

for n=1:1
%     bc_nodes=[bc_nodes-2*n; bc_nodes; bc_nodes+2*n];
end
bc=zeros(length(bc_nodes),2);
bc(:,1)=bc_nodes;


F=zeros(ndof,1);
F(last_nod)=last;

for n=1:1
    last_noder=[last_nod-2*n; last_nod+2*n];
    F(last_noder)=last/(3*n);
end

for n=1:1
    y_last=last_noder+1;
    F((y_last<last_nod)'*y_last)=n*last/10;
    F((y_last>last_nod)'*y_last)=-n*last/10;
end



% skapa K

K=zeros(ndof);
Ki=cell(nelm,1);

for el=1:nelm
    
    Ke=plani4e(Ex(el,:),Ey(el,:),ep,D);
    Ki{el}=Ke;
    indx=edof(el,2:end);
    K(indx,indx)=Ke+K(indx,indx);
    
end

A=ones(nelm,1);




%% SIMP

delta_0=1e-9; 
z=A;
K=K.*0;
for el=1:nelm
    indx=edof(el,2:end);
    K(indx,indx)=K(indx,indx)+(delta_0+(1-delta_0)*z(el)^p)*Ki{el};
end
Ks=sparse(K);
u=solveq(Ks,F,bc);



% g1=
g0save1=F'*u;
g1save1




%%
vM=ones(nelm,1);
dvM=ones(nelm,1);

P=[2 -1 0; -1 2 0; 0 0 6];


for el=1:nelm
    [es,~,eci]=plani4s(Ex(el,:),Ey(el,:),[1,t,1],D,ed(el,:));
    vm=sqrt(1/2*es*P*es');
    vM(el)=vm;
    dsig=es/z(el);
    dvM(el)=(1/vm)*P*es';
end


vM=(z.^0.5).*vM;
p_sig=8;
vMM=(sum(vM.^p_sig))^(1/p_sig);




%% Nummerisk derivata?? börja med gd0.
h=1e-3;
dg0n=zeros(nelm,1);


for el=1:nelm
    h_vec=zeros(nelm,1);
    h_vec(el)=h;
    zn=z.*(1+h_vec);
    indx=edof(el,2:end);
    Kn=K;
    Kn(indx,indx)=K(indx,indx)-(delta_0+(1-delta_0)*z(el)^p)*Ki{el}+(delta_0+(1-delta_0)*zn(el)^p)*Ki{el};
    Kns=sparse(Kn);
    un=solveq(Kns,F,bc);
    dg0n(el)=(F'*un-F'*u)/h;






end


%% Avståndsmatris.. (borde göra med Adjecency matris men nu har jag sparat för 2^5 iaf..)

M=zeros(nelm,nelm);

el_coord=[coord(enod(:,3),1)+coord(enod(:,2),1),coord(enod(:,5),2)+coord(enod(:,2),2)]/2;

for el=1:nelm
    for e=1:el
        M(el,e)=norm(el_coord(el,:)-el_coord(e,:));
    end
end

M=M+M';
granne=M<r;

% granne=sparse(granne);

w1=(1-M/r).*granne;
w2=sum((1-M/r).*granne)';



% fel derivata med MMA med svanberg



m = 1;
n = nelm;
z_t    = z;
xold1   = z;
xold2   = z;
xmin    = ones(nelm,1)*1e-6;
xmax    = ones(nelm,1);
low     = xmin;
upp     = xmax;
a0      = 1;
a       = zeros(m, 1);
c       = 1000*ones(m, 1);
d       = ones(m, 1);
outeriter = 0;
kkttol  = 1e-13;

x_mat=zeros(nelm,3);
x_mat(:,1)=z_t;
func_vals=zeros(1,3);
func_vals(1)=g0(u,F,z_t);

dg0_mat=zeros(nelm,1);

fig_nr=0;

for it=1:maxit
    outeriter=outeriter+1;
    dg0=zeros(nelm,1);

    for el=1:nelm
        ui=u(edof(el,2:end));
        t1=(3/(2*nelm)) * (sum(z_t)/nelm)^0.5 * F'*u;
        t2= p*(1-delta_0)*z_t(el)^(p-1)* (sum(z_t)/nelm)^1.5*ui'*Ki{el}*ui;
        dg0(el)=t1-t2;
    end
    
    dg0=w1./w2*dg0;
    
    dg0_mat(:,it)=dg0;

    g0val=g0(u,F,z_t)*1000;
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp]= ...
        mmasub(m,n,outeriter,z,xmin,xmax,xold1,xold2, ...
        g0val,dg0,g1val,dg1,low,upp,a0,a,c,d);
    xold2=xold1;
    xold1=z;
    z=xmma;



    w1=(1-M/r).*granne;
    w2=sum((1-M/r).*granne)';
    z_t=w1*z./w2;

    K=K.*0;
    for el=1:nelm
        indx=edof(el,2:end);
        K(indx,indx)=K(indx,indx)+(delta_0+(1-delta_0)*z_t(el)^p)*Ki{el};
    end


    Ks=sparse(K);
    u=solveq(Ks,F,bc);
    func_vals(it+1)=g0(u,F,z_t);
    if mod(outeriter,10)==0
        fig_nr=fig_nr+1;
        x_mat(:,fig_nr)=z_t;
    end

end


for fig=1:fig_nr
    figure()
    axis equal
    patch(Ex', Ey', x_mat(:,fig));

end

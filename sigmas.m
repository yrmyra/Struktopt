%% 

clc

L=0.01;
t=0.01;
E=4.5*10^9;
nu=0.32; 
V_box=12*L*2*L*t;
meshfac=2^0;
p=3;
last=-200;
maxit=3 ;




[coord, dof, enod, edof, Ex, Ey, bc] = designDomain(6*L, 2*L, L/meshfac);

nelm=length(edof);
nnod=length(coord);
ndof=nnod*2;

ep=[1 t 2];
D=E/(1-nu^2).*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];


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


% SIMP

delta_0=1e-9; 
z=A;
K=K.*0;
for el=1:nelm
    indx=edof(el,2:end);
    K(indx,indx)=K(indx,indx)+(delta_0+(1-delta_0)*z(el)^p)*Ki{el};
end
Ks=sparse(K);
u=solveq(Ks,F,bc);
ed=extract(edof,u);

g0(u,F,z);
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





%% vM derivata typ?





%%

m = 2;
n = nelm;
xval    = z;
xold1   = xval;
xold2   = xval;
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
x_mat(:,1)=xval;
func_vals=zeros(1,3);
func_vals(1)=F'*u;

ed_mat=cell(5,1);
fig_nr=0;

dg0_mat=zeros(nelm,1);
for it=1:maxit
    outeriter=outeriter+1;

    dg0=zeros(nelm,1);

    for el=1:nelm
        ui=u(edof(el,2:end));        
        dg0(el)=-p*(1-delta_0)*xval(el)^(p-1)*ui'*Ki{el}*ui;
    end

    dg0_mat(:,it)=dg0;

    g0val=F'*u;
    dg1=4*xval'/nelm;
    g1val=4*sum(xval)/nelm-1;
    g2val=vMM-1e5;
    dg2=sum(vM.^p_sig)^((1-p_sig)/p_sig)

    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp]= ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        g0val,dg0,g1val,dg1,low,upp,a0,a,c,d);
    xold2=xold1;
    xold1=xval;
    xval=xmma;

    for el=1:nelm
        grannar=[];
        noder=enod(el,2:5);
        for i=1:4
            [granne,~]=find(enod(:,2:5)==noder(i));
            grannar=unique([grannar;granne]);
        end
%         grannar(grannar==el)=[];
        w1=0;
        w2=0;
        for i=1:length(grannar)
            w=max(0,1-abs(z(el)-z(grannar(i))));
            w1=w1+w*z(grannar(i));
            w2=w2+w;
        end
        xval(el)=w1/w2;
        z(el)/xval(el);
        
    end

    K=K.*0;
    for el=1:nelm
        indx=edof(el,2:end);
        K(indx,indx)=K(indx,indx)+(delta_0+(1-delta_0)*xval(el)^p)*Ki{el};
    end


    Ks=sparse(K);
    u=solveq(Ks,F,bc);
    func_vals(it+1)=g0(u,F,xval);
    if mod(outeriter,10)==0
        fig_nr=fig_nr+1;
        ed=extract(edof,u);
        ed_mat{fig_nr}=ed;
        x_mat(:,fig_nr)=xval;
    end


end


for fig=1:fig_nr
    figure()
    patch(Ex', Ey', x_mat(:,fig));
end
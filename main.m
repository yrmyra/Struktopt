clc;clear;


% steg 1, hitta det mi ska minimera för grundformen
% - hitta compliance för grundfigur. Börja med ett F på något:))
% - vi ska minimera F'a*(V/Vbox)^(3/2)
% plan: tjocklek för alla emelent A(/x?). fan fattar inte hur man ska se
% elementen.. Fattar om man ser det om faktopr av densitet typ? Fattar om
% man ska optimera thickness...

%%
clc
clear

L=0.01;
t=0.01;
E=4.5*10^9;
nu=0.32; 
V_box=12*L*2*L*t;
meshfac=2^0;
plotpar=[ 1, 7, 1];




[coord, dof, enod, edof, Ex, Ey, bc] = designDomain(6*L, 2*L, L/meshfac);
% patch(Ex', Ey', 1)

% eldraw2(Ex,Ey,plotpar)


nelm=length(edof);
nnod=length(coord);
ndof=nnod*2;

A=ones(nelm,1)*L/10;
ep=[1 t 2];
D=E/(1-nu^2).*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];


% hitta noder med bc

last_nod=find(ismember(coord, [0.03, 0.02], 'rows'))*2;
% bc_nodes=find(ismember(coord, [0.005, 0;0.055, 0], 'rows'))*2;
bc_nodes=find(ismember(coord, [0.01, 0;0.05, 0], 'rows'))*2;


F=zeros(ndof,1);
F(last_nod)=-50;

% skapa K

K=zeros(ndof);
Ki=cell(nelm,1);

for el=1:nelm
    
    Ke=plani4e(Ex(el,:),Ey(el,:),ep,D);
    Ki{el}=Ke;
    K(edof(el,2:end),edof(el,2:end))=Ke+K(edof(el,2:end),edof(el,2:end));
    
end
% K=sparse(K);
bc=[bc_nodes(1) 0;bc_nodes(2) 0];

a=solveq(K,F,bc);

ed=extract(edof,a);
myeldisp2(Ex,Ey,ed,plotpar,10^3,A*10^3,1);


%% optimera comp och Vol?

% Vol just nu, se varje som konstant..? Testa optimera tjocklek till att
% börja med
%% SIMP först?

delta_0=1e-9; 
z=ones(nelm,1);
K=K.*0;
for el=1:nelm
    indx=edof(el,2:end);
    K(indx,indx)=K(indx,indx)+(delta_0+(1-delta_0)*z(el)^2)*Ki{el};
end

% K=sparse(K);
a=solveq(K,F,bc);

g0(a,F,z)

%% Nummerisk Diff


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

%%
dg0=zeros(nelm,1);
for el=1:nelm
    ai=a(edof(el,2:end));
    t1=(3/(2*nelm)) * (sum(z)/nelm)^0.5 * F'*a;
    t2= 2*(1-delta_0)*z(el)* (sum(z)/nelm)^1.5*ai'*Ki{el}*ai;
    dg0(el)=t1-t2;
%     dg0(el)=-ai'*Ki{el}*ai;
end

dg0./df

% (dg0/norm(dg0))./(df/norm(df))

%% Optimera map z
% testa exakt som innan trots att det typ är fel:))
% Tänker börja med CONLIN för gillar den mer

for it=1:3



end

















%% funktioner

function value=g0(a,F,z)
    value=F'*a*(sum(z)/length(z))^(3/2);
end








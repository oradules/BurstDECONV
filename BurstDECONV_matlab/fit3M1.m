function [res,resl,resh] = fit3M1( kmin, ksel,objmin,aa)
    

    %%%% compute 5 rates k1p,m k2p,m k3 from the 5 parameters %%%
    l1=kmin(1);
    l2=kmin(2);
    l3=kmin(3);
    A1=kmin(4);
    A2=kmin(5);
    A3=1-A1-A2;
    L1=l1+l2+l3;
    L2=l1.*l2+l1.*l3+l2.*l3;
    L3=l1.*l2.*l3;
    S1=A1.*l1+A2.*l2+A3.*l3;
    S2=A1.*l1.^2+A2.*l2.^2+A3.*l3.^2;
    S3=A1.*l1.^3+A2.*l2.^3+A3.*l3.^3;
    %%%% model M1
    k1p=-L3.*(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1p
    k2p=-(S2.^2-S1.*S3)./S1./(S1.^2-S2); %%% k2pk
    k3=-S1; %%%% k3
    k2m=(S1.^2-S2)./S1; %%% k2m
    k1m=-A1.*A2.*A3.*(l1-l2).^2.*(l1-l3).^2.*(l2-l3).^2.*S1./(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1m
    p1=k1m*k2m/(k1p*k2m+k1m*k2m+k1p*k2p);
    p2=k1p*k2m/(k1p*k2m+k1m*k2m+k1p*k2p);
    p3=k1p*k2p/(k1p*k2m+k1m*k2m+k1p*k2p);
    %%%% optimal
    res= [kmin(1),kmin(2),kmin(3),kmin(4),kmin(5),1-kmin(5),k1p,k1m,k2p,k2m,k3, p1,p2,p3,objmin,aa];
    %%%%% compute intervals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    l1=ksel(:,1);l1sel=l1;
    l2=ksel(:,2);l2sel=l2;
    l3=ksel(:,3);l3sel=l3;
    A1=ksel(:,4);A1sel=A1;
    A2=ksel(:,5);A2sel=A2;
    A3=1-A1-A2;A3sel=A3;
    L1=l1+l2+l3;
    L2=l1.*l2+l1.*l3+l2.*l3;
    L3=l1.*l2.*l3;
    S1=A1.*l1+A2.*l2+A3.*l3;S1sel=S1;
    S2=A1.*l1.^2+A2.*l2.^2+A3.*l3.^2;
    S3=A1.*l1.^3+A2.*l2.^3+A3.*l3.^3;
    %%%% model M1
    K1p=-L3.*(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1p
    K2p=-(S2.^2-S1.*S3)./S1./(S1.^2-S2); %%% k2p
    K3=-S1; %%%% k3
    K2m=(S1.^2-S2)./S1; %%% k2m
    K1m=-A1.*A2.*A3.*(l1-l2).^2.*(l1-l3).^2.*(l2-l3).^2.*S1./(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1m
    P1=K1m.*K2m./(K1p.*K2m+K1m.*K2m+K1p.*K2p);
    P2=K1p.*K2m./(K1p.*K2m+K1m.*K2m+K1p.*K2p);
    P3=K1p.*K2p./(K1p.*K2m+K1m.*K2m+K1p.*K2p);
    resl= [min(ksel(:,1)),min(ksel(:,2)),min(ksel(:,3)),min(ksel(:,4)),min(ksel(:,5)) ,min(1-ksel(:,5)),...
    min(K1p),min(K1m),min(K2p),min(K2m),min(K3),min(P1),min(P2),min(P3)];
    resh= [max(ksel(:,1)),max(ksel(:,2)),max(ksel(:,3)),max(ksel(:,4)),max(ksel(:,5)) ,max(1-ksel(:,5)),...
    max(K1p),max(K1m),max(K2p),max(K2m),max(K3),max(P1),max(P2),max(P3)];


end
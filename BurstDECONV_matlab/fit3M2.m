function [res,resl,resh] = fit3M2( kmin, ksel,objmin,aa)
    
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
    %%%%% Model M2
    kk3= -S1;
    kk1p = 1/2 * ( -L1+S2/S1 + sqrt((S1*L1-S2)^2-4*L3*S1)/S1 );
    kk2p = 1/2 * ( -L1+S2/S1 - sqrt((S1*L1-S2)^2-4*L3*S1)/S1 );
    kk1m = 1/2 * (S1-S2/S1 - (-S1^2*L1+S1*S2+S1*L2-L3+S2^2/S1-S3)/sqrt((S1*L1-S2)^2-4*L3*S1));
    kk2m = 1/2 * (S1-S2/S1 + (-S1^2*L1+S1*S2+S1*L2-L3+S2^2/S1-S3)/sqrt((S1*L1-S2)^2-4*L3*S1));
    pp1=kk1m*kk2p/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m);
    pp2=kk1p*kk2m/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m);
    pp3=kk1p*kk2p/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m);

    
    
    %%%% optimal
    res= [kmin(1),kmin(2),kmin(3),kmin(4),kmin(5),1-kmin(5),...
          kk1p,kk1m,kk2p,kk2m,kk3,pp1,pp2,pp3,objmin,aa];




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

    %%%% model M2
    KK1p = 1/2 * ( -L1+S2./S1 + sqrt((S1.*L1-S2).^2-4*L3.*S1)./S1 );
    KK2p = 1/2 * ( -L1+S2./S1 - sqrt((S1.*L1-S2).^2-4*L3.*S1)./S1 ); 
    KK1m = 1/2 * (S1-S2./S1 - (-S1.^2.*L1+S1.*S2+S1.*L2-L3+S2.^2./S1-S3)./sqrt((S1.*L1-S2).^2-4*L3.*S1)); 
    KK2m = 1/2 * (S1-S2./S1 + (-S1.^2.*L1+S1.*S2+S1.*L2-L3+S2.^2./S1-S3)./sqrt((S1.*L1-S2).^2-4*L3.*S1)); 
    KK3=-S1; %%%% k3
    PP1=KK1m.*KK2p./(KK1p.*KK2p+KK1m.*KK2p+KK1p.*KK2m);
    PP2=KK1p.*KK2m./(KK1p.*KK2p+KK1m.*KK2p+KK1p.*KK2m);
    PP3=KK1p.*KK2p./(KK1p.*KK2p+KK1m.*KK2p+KK1p.*KK2m);


    resl= [min(ksel(:,1)),min(ksel(:,2)),min(ksel(:,3)),min(ksel(:,4)),min(ksel(:,5)) ,min(1-ksel(:,5)),...
    min(KK1p),min(KK1m),min(KK2p),min(KK2m),min(KK3),min(PP1),min(PP2),min(PP3)];


    resh= [max(ksel(:,1)),max(ksel(:,2)),max(ksel(:,3)),max(ksel(:,4)),max(ksel(:,5)) ,max(1-ksel(:,5)),...
    max(KK1p),max(KK1m),max(KK2p),max(KK2m),max(KK3),max(PP1),max(PP2),max(PP3)];   

end
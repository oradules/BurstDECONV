function [res, resl, resh] = fit2(dirwrite,name,sd,dt,dtc,likelihood)
    %%%% define matrices that store parameters and objective functions
    store=[];
    fsz=16;
    if ~isempty(dt)

        for cens=0:1

            if cens
                [fs,xs,flo,fup]=ecdf([dt;dtc],'censoring',[zeros(size(dt));ones(size(dtc))]);
            else
                [fs,xs,flo,fup]=ecdf(dt);
            end
            %%%% fit distribution of spacings using combination of two exponentials
            xs=xs(1:end-1);
            fs=fs(1:end-1);
            flo=flo(1:end-1);
            fup=fup(1:end-1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            sN = sqrt(length(xs));

            %%%% try several cost functions
            if ~likelihood
                %%%% 1 : log scale only
                exp_fitness = @(k) (abs(log(k(3)*exp(k(1)*xs)+(1-k(3))*exp(k(2)*xs))-log(1-fs)))/sN; % k: parameters
                opts = optimoptions(@lsqnonlin,'TolFun', 1e-8,'MaxIter',1e3, ...
                'MaxFunEvals',1e6,'TolX', 1e-10);              
            else
                %%%% 3 : log likelihood, directly on spacings
                dt = [dt;dtc];
                options = optimoptions('fmincon','GradObj','on','TolFun', 1e-8,'TolX', 1e-10);
            end

            k00=[-0.01,-0.001,1];    
            amp = [log(100),log(100)]; 
            NbIterationinFit=100;    

            if likelihood
                lb=[-inf;-inf;0];ub=[0;0;1]; %%%% combined with fmincon
            end


            for mc = 1:NbIterationinFit

                %%%% Change k00
                factor=exp(amp.*(2*rand(1,2)-1)); 
                k0 = k00;
                k0(1:2)=k0(1,2).*factor;
                k0(3) =  2*rand(1,1)-1; %%% A1
                %%%% sort k0(1:2)
                k0(1:2)=sort(k0(1:2),'ascend');

                %%% impose constraints
                if ~( k0(1)*k0(3)+k0(2)*(1-k0(3)) < 0 )
                    while ~( k0(1)*k0(3)+k0(2)*(1-k0(3)) < 0 )
                        k0(3)=2*rand(1,1)-1;  %%% A1 
                    end
                end

                if ~likelihood        
                    %%%% Use the fcn lsqnonlin
                    [k, obj] = lsqnonlin(exp_fitness,k0,[],[],opts);
                    %%%% of use fmin con
                else
                   [k, obj] = fmincon(@maxlikelihood2withgradient,k0,[],[],[],[],lb,ub,[],options);
                end        
                %%%% sort k
                A=[k(3),1-k(3)]; %%% A values before sorting 
                [kk,IX]=sort(k(1:2),'ascend'); %%% sort lambdas
                k(1:2)=kk;
                A=A(IX);
                k(3)=A(1);

                store = [store;[k, obj, cens]];

                disp(mc)


            end    



        end % cens 


        %%%% select optimal and suboptimal obj < 2 objmin
        ind = find  (max(abs(imag(store(:,1:3))),[],2) < 1e-10);
        [objmin,indmin]=min(store(ind,4));
        imin=ind(indmin);
        overflow=1;
        ind=find(store(:,4) < (1+overflow)*objmin & max(abs(imag(store(:,1:3))),[],2) < 1e-10 );
        ksel=real(store(ind,1:3));
        kmin=real(store(imin,1:3) );
        censmin=store(imin,5);

        if censmin
            [fs,xs,flo,fup]=ecdf([dt;dtc],'censoring',[zeros(size(dt));ones(size(dtc))]);
        else
            [fs,xs,flo,fup]=ecdf([dt;dtc]);% no censoring
        end

        %%%%% plot survival function 
        h=figure(70)
        hold off
        semilogy(xs,1-fs,'or') %%% empirical function
        hold on
        semilogy(xs,1-flo,'--r') %%%% lower confidence 
        semilogy(xs,1-fup,'--r') %%%% upper confidence
        pred=kmin(3)*exp(kmin(1)*xs)+(1-kmin(3))*exp(kmin(2)*xs);
        semilogy(xs,pred,'k','linewidth',2) %%% predicted 2 exp
        axis([0, 250, 1e-6, 1])
        xlabel('Time [s]','fontsize',fsz)
        ylabel('Survival function','fontsize',fsz)


        %%%% compute 3 rates k1p,m k2 from the 3 parameters %%%
        S1 = kmin(3)*kmin(1)+(1-kmin(3))*kmin(2);
        k2 = -S1; 
        S2 = kmin(3)*(kmin(1))^2+(1-kmin(3))*(kmin(2))^2;
        S3 = kmin(3)*(kmin(1))^3+(1-kmin(3))*(kmin(2))^3;
        k1m = S1-S2/S1; 
        k1p = (S3*S1-S2^2)/S1/(S1^2-S2); 
        title(['k_1^-=',num2str(k1m,1),'k_1^+=',num2str(k1p,1),'k_2=',num2str(k2,1)])




        figfile=[dirwrite,'/Fit2_',name,'.pdf'];
        print(h,'-dpdf',figfile)



        %%%%%%% KS test

        thr= 20;
        ind = find(xs>thr); %%% perform the test on times larger than 20 s
        dist = max(abs( fs(ind)-1+pred(ind)));
        nsample=length(dt(dt>thr));
        N=10;
        c=sqrt(nsample)*dist;
        r=[1:N];
        aa=2*sum(exp(-2*c.^2*r.^2).*((-1).^(r-1)));

        h=figure(80)
        hold off
        semilogy(xs,fs,'kx')
        hold on
        semilogy(xs,1-pred,'ro')
        xlabel('Time [s]','fontsize',fsz)
        ylabel('CDF','fontsize',fsz)
        figfile=[dirwrite,'/Fit2_CDF_',name,'.pdf'];
        print(h,'-dpdf',figfile)
        %%%%%%%%%%%%%%%%%%%%%%%%%%






        %%%%% compute intervals
        S1 = ksel(:,3).*ksel(:,1)+(1-ksel(:,3)).*ksel(:,2);
        K2 = -S1; 
        S2 = ksel(:,3).*(ksel(:,1)).^2+(1-ksel(:,3)).*(ksel(:,2)).^2;
        S3 = ksel(:,3).*(ksel(:,1)).^3+(1-ksel(:,3)).*(ksel(:,2)).^3;
        K1m = S1-S2./S1; 
        K1p = (S3.*S1-S2.^2)./S1./(S1.^2-S2); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%         nn=size(DataExp);
        %%%% optimal
        %res= [k1p,k1m,k2,kmin(1),kmin(2),kmin(3),1-kmin(3),k1m/(k1m+k1p),k1p/(k1m+k1p),objmin,aa,nn(2),nn(1)];

        res= [kmin(1),kmin(2),kmin(3),1-kmin(3),k1p,k1m,k2,objmin,aa];

        P1=K1m./(K1m+K1p);
        P2=K1p./(K1m+K1p);

        resl= [min(ksel(:,1)),min(ksel(:,2)),min(ksel(:,3)),min(1-ksel(:,3)),min(K1p),min(K1m),min(K2)];
        resh= [max(ksel(:,1)),max(ksel(:,2)),max(ksel(:,3)),max(1-ksel(:,3)),max(K1p),max(K1m),max(K2)];
        
%         writecell({strrep(name,'result_','')},xlsfilename,'Sheet',1,'Range',['A',num2str(4*ifile-2)]); %%% filename
%         writematrix(res,xlsfilename,'Sheet',1,'Range',['B',num2str(4*ifile-2)]); %%% best result
%         writematrix(resl,xlsfilename,'Sheet',1,'Range',['B',num2str(4*ifile-1)]); %%% low
%         writematrix(resh,xlsfilename,'Sheet',1,'Range',['B',num2str(4*ifile)]); %%%  high


    else
        res = [0, 0, 0, 0, 0, 0, sd(2),sd(1)];
        resh= zeros(5);
        resl= zeros(5);
%         xlswrite(xlsfilename,{strrep(name,'result_','')},1,['A',num2str(4*ifile-2)]); %%% filename
%         nn=size(DataExp);
%         xlswrite(xlsfilename,res,1,['B',num2str(4*ifile-2)]); %%% best result

    end
end
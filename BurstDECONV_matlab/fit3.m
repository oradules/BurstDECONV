function [resM1, reslM1, reshM1,resM2, reslM2, reshM2] = fit3(dirwrite,name,sd,dt,dtc,likelihood)
    
    store=[];
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
                exp_fitness = @(k) (abs(log(k(4)*exp(k(1)*xs)+k(5)*exp(k(2)*xs)+(1-k(4)-k(5))*exp(k(3)*xs))-log(1-fs)))/sN; % k: parameters
                opts = optimoptions(@lsqnonlin,'TolFun', 1e-8,'MaxIter',1e3, ...
                'MaxFunEvals',1e6,'TolX', 1e-10);              
            else
                %%%% 3 : log likelihood, directly on spacings
                dt = [dt;dtc];
                options = optimoptions('fmincon','GradObj','on','TolFun', 1e-8,'TolX', 1e-10);
            end

            k00=[-0.1,-0.01,-0.001,0.25,0.25];    
            amp = [log(100),log(100),log(100),log(1),log(1)]; 
            NbIterationinFit=100;    

            if likelihood
               lb=[-inf;-inf;0];ub=[0;0;1]; %%%% combined with fmincon
            end


            for mc = 1:NbIterationinFit

                %%%% Change k00
                factor=exp(amp.*(2*rand(1,5)-1)); 

                k0 = k00.*factor;
                k0(4:5)=2*rand(1,2)-1; 

                %%%% sort k0(1:3)
                k0(1:3)=sort(k0(1:3),'ascend');

                if ~( sum(k0(4:5)) < 1 && k0(1)*k0(4)+k0(2)*k0(5)+k0(3)*(1-sum(k0(4:5))) < 0 )
                    while ~(sum(k0(4:5)) < 1 && k0(1)*k0(4)+k0(2)*k0(5)+k0(3)*(1-sum(k0(4:5))) < 0)
                        k0(4:5)=2*rand(1,2)-1;  %%% A1,A2 values
                    end
                end    

                if ~likelihood        
                    %%%% Use the fcn lsqnonlin
                    [k, obj] = lsqnonlin(exp_fitness,k0,[],[],opts);
                else
                    [k, obj] = fmincon(@maxlikelihood3withgradient,k0,[],[],[],[],lb,ub,[],options);
                end        
                %%%% sort k
                A=[k(4),k(5),1-k(4)-k(5)]; %%% A values before sorting 
                [kk,IX]=sort(k(1:3),'ascend'); %%% sort lambdas
                k(1:3)=kk;
                A=A(IX);
                k(4:5)=A(1:2);

                store = [store;[k, obj, cens]];

                disp([num2str(mc),'/',num2str(NbIterationinFit)])


            end    



        end % cens 


   %%%%%%%%%%%%%%%%%%%%%%%%% parametric fit %%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fsz=16;
    %%%% select optimal and suboptimal obj < 2 objmin
    ind = find  (max(abs(imag(store(:,1:3))),[],2) < 1e-10); % takes the part where we have the imaginary part of lambda i almost 0
    [objmin,indmin]=min(store(ind,end-1)); %minimum of the objectives for all the previous lambdas
    imin=ind(indmin); % finding the positions of the lambda where we have the lowest objective
    overflow=1;
    ind=find(store(:,end-1) < (1+overflow)*objmin & max(abs(imag(store(:,1:3))),[],2) < 1e-10 ); % taking suboptimal objectives such that they are <2*optimal objective and they satisfy the fact that the imaginary part of lambda is almost 0
    ksel=real(store(ind,1:5)); % taking an array of lambda i's and Ai's where we have the the objective function less than <2*times the minimum objective function
    kmin=real(store(imin,1:5) ); %taking the lambda i's and A i's that are the fittest
    censmin=store(imin,end);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if censmin
        [fs,xs,flo,fup]=ecdf([dt;dtc],'censoring',[zeros(size(dt));ones(size(dtc))]);
    else
        [fs,xs,flo,fup]=ecdf([dt;dtc]);% no censoring
    end
    %%%%% plot survival function 
    h=figure(70);
    hold off
    semilogy(xs,1-fs,'or') %%% empirical function
    hold on
    semilogy(xs,1-flo,'--r') %%%% lower confidence 
    semilogy(xs,1-fup,'--r') %%%% upper confidence
    pred=kmin(4)*exp(kmin(1)*xs)+kmin(5)*exp(kmin(2)*xs)+(1-kmin(4)-kmin(5))*exp(kmin(3)*xs);
    semilogy(xs,pred,'k','linewidth',2) %%% predicted 2 exp
    axis([0, 250, 1e-6, 1])
    xlabel('Time [s]','fontsize',fsz)
    ylabel('Survival function','fontsize',fsz)
    figfile=[dirwrite,'/Fit3_',name,'.pdf'];
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
    h=figure(80);
    hold off
    semilogy(xs,fs,'kx')
    hold on
    semilogy(xs,1-pred,'ro')
    xlabel('Time [s]','fontsize',fsz)
    ylabel('CDF','fontsize',fsz)
    figfile=[dirwrite,'/Fit3_CDF_',name,'.pdf'];
    print(h,'-dpdf',figfile)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        [resM1,reslM1,reshM1]=fit3M1( kmin, ksel,objmin,aa);
        [resM2,reslM2,reshM2]=fit3M2( kmin, ksel,objmin,aa);

    else
        resM1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, sd(2), sd(1)];
        reslM1= zeros(8,1);
        reshM1= zeros(8,1);
        resM2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, sd(2), sd(1)];
        reslM2= zeros(8,1);
        reshM2= zeros(8,1);
        
    end
    
end
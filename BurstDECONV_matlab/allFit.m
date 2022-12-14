%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier
%%%%% Copyright : This is published under 3-clause BSD
%%%%% last change november 2021
%%%%% infers multiple models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


global FreqEchSimu FreqEchImg DureeAnalysee TaillePreMarq ...
            TailleSeqMarq TaillePostMarq  Polym_speed frame_num num_possible_poly EspaceInterPolyMin ...
            DureeSimu Intensity_for_1_Polym;
global dt;
tic

%%%% load parameters
[ Polym_speed, ~,TaillePreMarq,TailleSeqMarq,TaillePostMarq,EspaceInterPolyMin,FrameLen,Intensity_for_1_Polym] = parameters();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FreqEchImg = (1/FrameLen); % image per second data time sampling   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed); % how many interval(possible poly start position) in 1s


fsz=16;lw=2;
sz=10;


likelihood=0; %%%%% if this is true use max likelihood for parametric survival function fit

if likelihood
    DataFilePath0 = 'output/FitResults_likelihood'; %%%%% where to write fit results
else
    DataFilePath0 = 'output/FitResults'; %%%%% where to write fit results    
end
mkdir(DataFilePath0);
DataFilePath='output/matfiles/'; %%%%% where are the deconvolution results

ifile=1;

%%% list of movies to process
lnames={'dataD32Sx'};
%%% or process all result files as the same phenotype
data_list=struct2cell(dir(fullfile(DataFilePath)));
file_name_list = data_list(1,3:end);
lnames=[];
for i=1:length(file_name_list)
    if strfind(file_name_list{i},'result')
       lnames = [lnames,{file_name_list{i}}];
    end
end    
%%%% name of the phenotype
name=strrep(strrep(lnames{1},'.mat',''),'result_',''); %%%%%%%%%%
dirwrite=[DataFilePath0,'/',name,'_result'];
mkdir(dirwrite);
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% creating xls result files for all types of models      
%%% 2 states
xlsfilename2 = [dirwrite,'/fit2_results.xlsx'];
xlswrite(xlsfilename2,{'Data','lambda1','lambda2','A1','A2','k1p','k1m','k2','Obj','KS test'},1,'A1');
%%% 3 states M1
xlsfilename3M1 = [dirwrite,'/fit3_M1_results.xlsx'];
xlswrite(xlsfilename3M1,{'Data','lambda1','lambda2','lambda3','A1','A2','A3',...
    'k1p','k1m','k2p','k2m','k3','p1','p2','p3',...
    'Obj','KS test'},1,'A1')
%%% 3 states M1
xlsfilename3M2 = [dirwrite,'/fit3_M2_results.xlsx'];
xlswrite(xlsfilename3M2,{'Data','lambda1','lambda2','lambda3','A1','A2','A3',...
    'k1p','k1m','k2p','k2m','k3','p1','p2','p3',...
    'Obj','KS test'},1,'A1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(lnames) > 1
%%%%% lump files, Procustes method 
dataExp=[];
dataPred=[];
posPred=[];
tmax=[];
%%%% first compute max dimension
nmax=1;nmaxpos=1;
for iii=1:length(lnames)
    fname = [DataFilePath,lnames{iii}];
    load(fname);
    n2=size(DataExp);
    n3=size(PosPred);
    if n2(1) > nmax
        nmax=n2(1);
    end
    if n3(1) > nmaxpos
        nmaxpos=n3(1);
    end
end

for iii=1:length(lnames)
    fname = [DataFilePath,lnames{iii}];
    load(fname);
    n2=size(DataExp);
    n3=size(PosPred);    
%%%%%%%%%%%%%% padding 
    DataExp=[DataExp;zeros(nmax-n2(1),n2(2))];
    DataPred=[DataPred;zeros(nmax-n2(1),n2(2))]; 
    PosPred=[PosPred;zeros(nmaxpos-n3(1),n3(2))];
%%%%%%%%%%%%% concatenating
    dataExp=[dataExp,DataExp]; 
    dataPred=[dataPred,DataPred];
    posPred=[posPred,PosPred];    
    tmax=[tmax;(n2(1)-1)/FreqEchImg*ones(n2(2),1)]; %%%% movie length            
end
%%%%% padded and concatenated data%%%%
    DataExp=dataExp;
    DataPred=dataPred;
    PosPred=posPred;
else
    fname = [DataFilePath,lnames{1}];
    load(fname);   
    n2=size(DataExp);
    n3=size(PosPred);
    tmax= (n2(1)-1)/FreqEchImg*ones(n2(2),1); %%%% end of the movie for all cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% if only one file per phenotype this is just the content of the file
%%%% with no padding and concatenation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%% if length(lnames) > 1
    
[dt,dtc,T0] = dtcomp(DataExp,DataPred,PosPred,tmax,FrameLen,EspaceInterPolyMin,Polym_speed,TaillePostMarq,TaillePreMarq,TailleSeqMarq);
n=size(DataExp);

%%%%%%% show and print data and deconvolution fit %%%%%%%%
%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iii=1:n(2) %%% for all cells    
    ifig=floor((iii-1)/36)+1;
    h=figure(ifig);
    hold off
    subplot(6,6,mod((iii-1),36)+1)
    hold off
    plot((0:(n(1)-1))/FreqEchImg,DataExp(:,iii),'k','linewidth',0.1)
    hold on
    prediction=sumSignal1_par(find((PosPred(:,iii)==1))',FreqEchSimu, FreqEchImg, TaillePreMarq, ...
            TailleSeqMarq, TaillePostMarq,  Polym_speed, n(1), Intensity_for_1_Polym);
    plot((0:(n(1)-1))/FreqEchImg,prediction,'r','linewidth',0.1)
    %pos= [T0(iii),0,tmax(iii)-T0(iii),100];
    pos =[T0(iii),0; tmax(iii),0;  tmax(iii),100; T0(iii),100];
    %rectangle('Position', pos, 'FaceColor', [0 0 1 0.1], 'EdgeColor',[0 0 1 0.1]); 
    patch('vertices',pos,'faces',[1,2,3,4],'FaceColor','b','FaceAlpha',0.1,'LineStyle','none'); %%% for older Matlab versions

    axis([0,tmax(iii),0,100])
    title(num2str(iii))  
    if mod(iii,36)==0 || iii == n(2)
        figfile = [dirwrite,'/','fig',num2str(ifig),'.pdf']
    print(h,'-dpdf',figfile) 
    end
end %%% for iii

%%%%%%%%%%%%% save parameters results for 2 state model
[res, resl, resh] = fit2(dirwrite,name,n,dt,dtc,likelihood);

xlswrite(xlsfilename2,{strrep(name,'result_','')},1,['A',num2str(4*ifile-2)]); %%% filename
xlswrite(xlsfilename2,res,1,['B',num2str(4*ifile-2)]); %%% best result
xlswrite(xlsfilename2,resl,1,['B',num2str(4*ifile-1)]); %%% low
xlswrite(xlsfilename2,resh,1,['B',num2str(4*ifile)]); %%%  high

%%%%%%%%%%%% save parameters results for 3 state model            
[resM1, reslM1, reshM1,resM2, reslM2, reshM2]=fit3(dirwrite,name,n,dt,dtc,likelihood);


%%%%% Model M1
xlswrite(xlsfilename3M1,{strrep(name,'result_','')},1,['A',num2str(4*ifile-2)]); %%% filename
xlswrite(xlsfilename3M1,resM1,1,['B',num2str(4*ifile-2)]); %%% best result
xlswrite(xlsfilename3M1,reslM1,1,['B',num2str(4*ifile-1)]); %%% low
xlswrite(xlsfilename3M1,reshM1,1,['B',num2str(4*ifile)]); %%%  high
%%%%% Model M2
xlswrite(xlsfilename3M2,{strrep(name,'result_','')},1,['A',num2str(4*ifile-2)]); %%% filename
xlswrite(xlsfilename3M2,resM2,1,['B',num2str(4*ifile-2)]); %%% best result
xlswrite(xlsfilename3M2,reslM2,1,['B',num2str(4*ifile-1)]); %%% low
xlswrite(xlsfilename3M2,reshM2,1,['B',num2str(4*ifile)]); %%%  high



toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% Copyright : This is published under 3-clause BSD
%%%%% Last change March 2021
%%%%% reads data from excel files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%% all the excel data files should be placed here
DataFilePath = 'data/';
%%% all output files will be placed here
OutputPath = 'output/';
mkdir(OutputPath);

%%% all the generated mat files will be placed here
MatFilePath = [OutputPath,'matfiles/'];
mkdir(MatFilePath);
%%%% plots of traces will be placed here
DataFilePath0 = [OutputPath,'images/'];
mkdir(DataFilePath0);

data_list=struct2cell(dir(fullfile(DataFilePath)));
file_name_list = data_list(1,3:end);
nexp=length(file_name_list); 


for data_i=1:nexp
%%%%%% where the data is %%%%%%%%%%%%%%%%%%%
DataFileName = [DataFilePath,file_name_list{data_i}];
%%%%the fields below should cover data and excess should be empty
[num,Spot_names]=xlsread(DataFileName,1,'E1:ON1');
DataExp = xlsread(DataFileName,1,'E3:ON700'); %%%% up to 400 spots

n=size(DataExp);
fn=strrep(strrep(strrep(file_name_list{data_i},'_',''),'.xlsx',''),'.xlsx','');

sd= size(DataExp);
Frames=0:sd(1)-1;
Samples=1:sd(1);
fname = [MatFilePath,'data_',fn,'.mat'];
save(fname,'DataExp')


for iii=1:n(2)    
    ifig=floor((iii-1)/36)+1;
    h=figure(ifig);
    hold off
    subplot(6,6,mod((iii-1),36)+1)
    plot(1:n(1),DataExp(:,iii),'k','linewidth',0.1)
    %axis([1 500 0 100])
    title(num2str(iii))  
    title(strrep(Spot_names{iii},'_','-'))
    if mod(iii,36)==0 || iii == n(2)
        figfile = [DataFilePath0,fn,'_','fig',num2str(ifig),'.pdf']
    print(h,'-dpdf',figfile) 
    end


end %%% for iii



end %%% data_i


name_list=file_name_list;
for i=1:length(file_name_list)
    name_list{i}=[num2str(i,3),': ',name_list{i}];
end
%%%%% write text file
fname = [OutputPath,'description.txt']
tfile = fopen(fname,'w')

fprintf( tfile, '%s \n',name_list{:})
fclose(tfile)
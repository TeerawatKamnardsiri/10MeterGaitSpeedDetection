%Teerawat Kamnardsiri
%Colledge of Arts, Media and Technology, Chiang Mai University
%Loading data of Approach Run
%ComparingGS and Mocap split of walking speed.
%December, 2022.
clc;
close all;
clear;

%Output results data
Outputpath = '/ResultsComparingTwoSystems/';

%Video Capture data path
Videopath = '/ResultsVideo/';
VideoName='S1S5N2Speed';

%Motion Capture data path
MoCappath = '/ResultsMocap/';
MoCapName='S1S5N2_CleanedSpeed';

%************************************************************************
% Gait Speed detection System
% Load raw data from exel file%Gait speed detection

Video_File=fullfile([VideoName,'.xlsx']);

try
VideoData = fullfile([Videopath,Video_File]);
rawVideoData = xlsread(VideoData);
% Defile point 1
catch
    fprintf('File of the data not found, the directory should be changed!');
    break;
end

% Defile point of Video data type
GaitDataX=rawVideoData(:,1);
GaitDataV=rawVideoData(:,2);

% Setting filter parameters 2th order, cut of sample 30 Hz, 
% cutoff frequency 30/2 = 15 as 6/15 = 0.4 of Nyquist frequency.
% cutoff frequency 60/2 = 30 as 6/30 = 0.2 of Nyquist frequency.
% cutoff frequency 120/2 = 60 as 6/60 = 0.1 of Nyquist frequency.

[b,a]=butter(2,0.2);
ButterGaitDataV = filter(b,a,GaitDataV);
%ButterGaitDataX = filter(b,a,GaitDataX);

ButterGaitDataV = smooth(GaitDataV);
%GaitDataX = smooth(GaitDataX);

%plot(GaitDataV,GaitDataX)

%Setting index of each period 
indGaitsix=1;
indGaitfive=1;
indGaitfour=1;
indGaitthree=1;
indGaittwo=1;
indGaitone=1;

    
for Gi=1:size(ButterGaitDataV)
    
    % period 6-5 m
    if  GaitDataX(Gi)<6.00 && GaitDataX(Gi)>=5.00 
        NewGaitV_six(indGaitsix)=ButterGaitDataV(Gi);
        NewGaitX_six(indGaitsix)=GaitDataX(Gi);
        indGaitsix=indGaitsix+1;
        
    % period 5-4 m
    elseif  GaitDataX(Gi)<5.00 && GaitDataX(Gi)>=4.00 
        NewGaitV_five(indGaitfive)=ButterGaitDataV(Gi);
        NewGaitX_five(indGaitfive)=GaitDataX(Gi);
        indGaitfive=indGaitfive+1;
    
    % period 4-3 m    
    elseif  GaitDataX(Gi)<4.00 && GaitDataX(Gi)>=3.00 
        NewGaitV_four(indGaitfour)=ButterGaitDataV(Gi);
        NewGaitX_four(indGaitfour)=GaitDataX(Gi);
        indGaitfour=indGaitfour+1;
   
        % period 3-2 m    
    elseif  GaitDataX(Gi)<3.00 && GaitDataX(Gi)>=2.00 
        NewGaitV_three(indGaitthree)=ButterGaitDataV(Gi);
        NewGaitX_three(indGaitthree)=GaitDataX(Gi);
        indGaitthree=indGaitthree+1;

        % period 2-1 m    
    elseif  GaitDataX(Gi)<2.00 && GaitDataX(Gi)>=1.00 
        NewGaitV_two(indGaittwo)=ButterGaitDataV(Gi);
        NewGaitX_two(indGaittwo)=GaitDataX(Gi);
        indGaittwo=indGaittwo+1;
        
        % period 1-0 m    
    elseif  GaitDataX(Gi)<1.00 && GaitDataX(Gi)>=0.00 
        NewGaitV_one(indGaitone)=ButterGaitDataV(Gi);
        NewGaitX_one(indGaitone)=GaitDataX(Gi);
        indGaitone=indGaitone+1;    
    end
end

    % Mean, SD, COV of period 6-5 m
    Gait_Mean6= mean(NewGaitV_six);
    %Gait_SD6= std(NewGaitV_six);
    
    % Mean, SD, COV of period 5-4 m
    Gait_Mean5= mean(NewGaitV_five);
    %Gait_SD5= std(NewGaitV_five);
    

    % Mean, SD, COV of period 4-3 m
    Gait_Mean4= mean(NewGaitV_four);
    %Gait_SD4= std(NewGaitV_four);
    
    
    % Mean, SD, COV of period 3-2 m
    Gait_Mean3= mean(NewGaitV_three);
    %Gait_SD3= std(NewGaitV_three);
    
    
    % Mean, SD, COV of period 2-1 m
    Gait_Mean2= mean(NewGaitV_two);
    %Gait_SD2= std(NewGaitV_two);
    
    % Mean, SD, COV of period 1-0 m
    Gait_Mean1= mean(NewGaitV_one);
    %Gait_SD1= std(NewGaitV_one);


%************************************************************************
% 3D Motion Capture System
% Load raw data from exel file
% filename='T01_FS01cropV_pre';

MoCapData = fullfile([MoCappath,MoCapName,'.xlsx']);
rawMoCapData = xlsread(MoCapData);

% Defile point 1
MoCapTime=rawMoCapData(:,1);
MoCapX=rawMoCapData(:,2);
MoCapV=rawMoCapData(:,3);

% Setting filter parameters 2th order, cut of sample 30 Hz, 
% cutoff frequency 30/2 = 15 as 6/15 = 0.4 of Nyquist frequency.
% cutoff frequency 60/2 = 30 as 6/30 = 0.2 of Nyquist frequency.
% cutoff frequency 120/2 = 60 as 6/60 = 0.1 of Nyquist frequency.


%MoCap Data filtering
[d,c]=butter(2,0.1);
ButterMoCapV = filter(d,c,MoCapV);
ButterMoCapX = filter(d,c,MoCapX);

%ButterMoCapV = (MoCapV);

%Setting index of each period 
indMocapsix=1;
indMocapfive=1;
indMocapfour=1;
indMocapthree=1;
indMocaptwo=1;
indMocapone=1;

% period 5-4 m
for Mi=1:size(ButterMoCapV)
    % period 6-5 m
    if MoCapX(Mi)<6.00 && MoCapX(Mi)>=5.00 
        NewMoCapV_six(indMocapsix)=ButterMoCapV(Mi);
        NewMoCapX_six(indMocapsix)=MoCapX(Mi);
        indMocapsix=indMocapsix+1;
    % period 5-4 m
    elseif MoCapX(Mi)<5.00 && MoCapX(Mi)>=4.00 
        NewMoCapV_five(indMocapfive)=ButterMoCapV(Mi);
        NewMoCapX_five(indMocapfive)=MoCapX(Mi);
        indMocapfive=indMocapfive+1;
   % period 4-3 m     
        elseif  MoCapX(Mi)<4.00 && MoCapX(Mi)>=3.00 
            NewMoCapV_four(indMocapfour)=ButterMoCapV(Mi);
            NewMoCapX_four(indMocapfour)=MoCapX(Mi);
            indMocapfour=indMocapfour+1;
   % period 3-2 m          
        elseif  MoCapX(Mi)<3.00 && MoCapX(Mi)>=2.00 
            NewMoCapV_three(indMocapthree)=ButterMoCapV(Mi);
            NewMoCapX_three(indMocapthree)=MoCapX(Mi);
            indMocapthree=indMocapthree+1;
   % period 2-1 m          
        elseif  MoCapX(Mi)<2.00 && MoCapX(Mi)>=1.00 
            NewMoCapV_two(indMocaptwo)=ButterMoCapV(Mi);
            NewMoCapX_two(indMocaptwo)=MoCapX(Mi);
            indMocaptwo=indMocaptwo+1;          
      % period 1-0 m          
        elseif  MoCapX(Mi)<1.00 && MoCapX(Mi)>=0.00 
            NewMoCapV_one(indMocapone)=ButterMoCapV(Mi);
            NewMoCapX_one(indMocapone)=MoCapX(Mi);
            indMocapone=indMocapone+1;    
    end
end
    % Mean, SD, COV of period 6-5 m
    Mocap_Mean6= mean(NewMoCapV_six);
    %Mocap_SD6= std(NewMoCapV_six);
    %Mocap_COV6= (Mocap_SD6/Mocap_Mean6)*100;
    
    % Mean, SD, COV of period 5-4 m
    Mocap_Mean5= mean(NewMoCapV_five);
    %Mocap_SD5= std(NewMoCapV_five);
    %Mocap_COV5= (Mocap_SD5/Mocap_Mean5)*100;
    
    % Mean, SD, COV of period 4-3 m
    Mocap_Mean4= mean(NewMoCapV_four);
    %Mocap_SD4= std(NewMoCapV_four);
    %Mocap_COV4= (Mocap_SD4/Mocap_Mean4)*100;
    
    % Mean, SD, COV of period 3-2 m
    Mocap_Mean3= mean(NewMoCapV_three);
    %Mocap_SD3= std(NewMoCapV_three);
    %Mocap_COV3= (Mocap_SD3/Mocap_Mean3)*100;   
    
    % Mean, SD, COV of period 2-1 m
    Mocap_Mean2= mean(NewMoCapV_two);
    %Mocap_SD2= std(NewMoCapV_two);
    %Mocap_COV2= (Mocap_SD2/Mocap_Mean2)*100;  

    % Mean, SD, COV of period 1-0 m
    Mocap_Mean1= mean(NewMoCapV_one);
    %Mocap_SD1= std(NewMoCapV_one);
    %Mocap_COV1= (Mocap_SD1/Mocap_Mean1)*100;  
    
        
    
fprintf('            \tMocap \tGait Speed\n');
fprintf('___________________________\n');
fprintf('Mean 0-1    \t%0.3f \t%0.3f \n',Mocap_Mean6,Gait_Mean6);
%fprintf('SD6      \t%0.3f \t%0.3f \n',Mocap_SD6,Gait_SD6);
%fprintf('COV6     \t%0.3f \t%0.3f\n',Mocap_COV6,Gait_COV6);
fprintf('___________________________\n');
fprintf('Mean 1-2    \t%0.3f \t%0.3f \n',Mocap_Mean5,Gait_Mean5);
%fprintf('SD5      \t%0.3f \t%0.3f \n',Mocap_SD5,Gait_SD5);
%fprintf('COV5     \t%0.3f \t%0.3f\n',Mocap_COV5,Gait_COV5);
fprintf('___________________________\n');
fprintf('Mean 2-3    \t%0.3f \t%0.3f \n',Mocap_Mean4,Gait_Mean4);
%fprintf('SD4      \t%0.3f \t%0.3f \n',Mocap_SD4,Gait_SD4);
%fprintf('COV4     \t%0.3f \t%0.3f\n',Mocap_COV4,Gait_COV4);
fprintf('___________________________\n');
fprintf('Mean 3-4    \t%0.3f \t%0.3f \n',Mocap_Mean3,Gait_Mean3);
%fprintf('SD3      \t%0.3f \t%0.3f \n',Mocap_SD3,Gait_SD3);
%fprintf('COV3     \t%0.3f \t%0.3f\n',Mocap_COV3,Gait_COV3);
fprintf('___________________________\n');
fprintf('Mean 4-5    \t%0.3f \t%0.3f \n',Mocap_Mean2,Gait_Mean2);
%fprintf('SD2      \t%0.3f \t%0.3f \n',Mocap_SD2,Gait_SD2);
%fprintf('COV2     \t%0.3f \t%0.3f\n',Mocap_COV2,Gait_COV2);
fprintf('___________________________\n');
fprintf('Mean 5-6    \t%0.3f \t%0.3f \n',Mocap_Mean1,Gait_Mean1);
%fprintf('SD1      \t%0.3f \t%0.3f \n',Mocap_SD1,Gait_SD1);
%fprintf('COV1     \t%0.3f \t%0.3f\n',Mocap_COV1,Gait_COV1);



%subplot(2,1,1)
%figure1 = figure;
%hold all;

%plot(MoCapX, ButterMoCapV);
%plot(GaitDataX, ButterGaitDataV);

% legend('Mo-Cap','Gait Speed','Mo-Cap (Vmax)','Gait Speed (Vmax)','Mo-Cap (Vmin)','Gait Speed (Vmin)');
%legend('Mo-Cap','Gait Speed');
%xlabel('X Data (m)');
%ylabel('Velocity Data (m/s)');


%Saving the relationship between Position and Speed
%saveas(figure1,(fullfile([datapath,filename,'Gait_and_MoCap_plot.jpg'])));
 


% %Saving data to the excel file
% %A = {'MoCap_Vmax','GS_Vmax','MoCap_Vmin','GS_Vmin','DMoCap_Vmax','DGS_Vmax','DMoCap_Vmin','DGS_Vmin'; VmaxpksT, VmaxpksG, -VminpksT, -VminpksG, 6-MoCapX(VmaxlocsT), 6-GaitDataX(VmaxlocsG),6-MoCapX(VminlocsT),6-GaitDataX(VminlocsG)};
 A = {'MoCap(0-1)','GS_Mean(0-1)','MoCap(1-2)','GS_Mean(1-2)','MoCap_Mean(2-3)','GS_Mean(2-3)','MoCap_Mean(3-4)','GS_Mean(3-4)','MoCap_Mean(4-5)','GS_Mean(4-5)','MoCap_Mean(5-6)','GS_Mean(5-6)'; Mocap_Mean6, Gait_Mean6, Mocap_Mean5, Gait_Mean5, Mocap_Mean4, Gait_Mean4, Mocap_Mean3, Gait_Mean3, Mocap_Mean2, Gait_Mean2,Mocap_Mean1, Gait_Mean1 };
 sheet = 1;
 xlRange = 'A1';
 xlswrite(fullfile([Outputpath,VideoName,'Mean_Results.xlsx']),A,sheet,xlRange);
 
fprintf('***************************************************************\n');
fprintf('Data has been successfully saved.');

% promptMessage = sprintf('Data has been successfully saved.');
% button = questdlg(promptMessage, 'OK', 'OK', 'Cancel', 'OK');
% if strcmpi(button, 'Cancel')
%   return; % Or break or continue
% end

close All;

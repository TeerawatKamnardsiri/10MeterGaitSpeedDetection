%Teerawat Kamnardsiri
%Colledge of Arts, Media and Technology, Chiang Mai University
%Gait Speed Detection System
%December, 2022.
clc;
close all;
clear;

datapath = 'S1/';
Outputpath = 'S1/';
BGFileName='BG';

FileName='S1S1F1';

BG_File=fullfile([BGFileName,'.MOV']);
FG_File=fullfile([FileName,'.MOV']);

%Distance of the camera around 18-20 m., Sony5100 with Lens 20 mm.
PIXAREA=700;

%Setting Threshold
THRESHOLDSET=30;

%Calibration distance setup
MARKER_CALIBRATION=6;

%Set Motion Blur 
MOTIONBLUR=30;

%Foreground starting frame
StartCount=120;

%Cropping Area of a image 
Crop = importdata('position1.txt');
P1 = Crop(:,1);
P2 = Crop(:,2);
P3 = Crop(:,3);
P4 = Crop(:,4);

%*****************************************************************************************************************
%Read Background Image
BG= VideoReader(sprintf('%s%s',datapath,BG_File));
%Read Background Image
Background = read(BG,2);
Background = imresize(Background, 0.5);

I=Background;
I1 = rgb2gray(I) ;
[y,x] = find(I1) ;
R = I(:,:,1) ; G = I(:,:,2) ; B = I(:,:,3) ;
x0 = min(x) ; x1 = max(x) ; 
y0 = min(y) ; y1 = P2 ;

[X1,Y1] = meshgrid(x0:x1,y0:y1) ;
x2 = min(x) ; x3 = max(x); 
y2 = (P2+P4) ; y3 = max(y);

[X2,Y2] = meshgrid(x2:x3,y2:y3) ;

idx1 = sub2ind(size(I1),Y1(:),X1(:)) ;
idx2 = sub2ind(size(I1),Y2(:),X2(:)) ;
I2 = I ;

%C = [mean(R(idx)) mean(G(idx)) mean(B(idx))] ;
C = [R(255),G(255),B(255)];
for i = 1:3
    T2 = I2(:,:,i) ;
    T2(idx1) = C(i) ;
    T2(idx2) = C(i) ;
    I2(:,:,i) = T2 ;
end




%Size of Background2
Background2W = BG.Width*0.5;
Background2H = BG.Height*.05;

Background_BW = rgb2gray(Background);
Background_BW = imgaussfilt(Background_BW,2);
Background_BW = imadjust(Background_BW,[0.00; 0.4],[0.1; 1.0], 0.35); 
    
%Calibration the system
box_color = {'white'};
Title = sprintf(' 6-Meter Gait Speed Detection System \n \n Calibration the system, please perform 2 steps:\n 1: click mouse at the Marker M1 and then\n 2: click mouse at the Marker M2 (10 m. from M1) for calibration of the system.');

Background2 = insertText(Background, [Background2W/2 Background2H/100], Title, 'FontSize',18,'BoxColor',...
    box_color,'BoxOpacity',0.5,'TextColor','black','AnchorPoint','CenterTop');

imshow(Background2);

%imshow(Background_BW);

%Loading video
Gait= VideoReader(sprintf('%s%s',datapath,FG_File));
GaitW= Gait.Width*0.5;
GaitH= Gait.Height*0.5;

Foreground=read(Gait,1);
Foreground = imresize(Foreground, 0.5);
Foreground_BW = rgb2gray(Foreground);
Foreground_BW = imadjust(Foreground_BW,[0.00; 0.4],[0.1; 1.0], 0.35);   

%Getting  area positions for processing

position = dlmread('position1.txt');

%Getting the positions of the start point and the end point
[PositionX,PositionY]=ginput(2);

%Getting the distance of the Approach-Run.

EndPoint = PositionX(1);
TenMeters = PositionX(2)- PositionX(1);

StartPoint = EndPoint+TenMeters;

OneMeterW=(TenMeters/MARKER_CALIBRATION); %pixels per meter Width
OnePixelW=(MARKER_CALIBRATION/TenMeters); %meter per pixel Width

MaxDistance=GaitW/(TenMeters/MARKER_CALIBRATION);
%Distance at the Start point to the End pointconvert to pixel format in the image.
DistanceStartToEnd = StartPoint-EndPoint;
RealTotalDistance = DistanceStartToEnd * OnePixelW;

%***************************************************************************************************************
%Loadding Foreground video
DurationTime = 4;% The value should be (4,8,10,12,16,20,24,32)
StepTime = 4;
%FrameRate = 60;
FrameRate=Gait.FrameRate;
Velocity = 0;
PositionBody=0;

%Read and Starting at specific time
NumberofFrame= ceil(Gait.FrameRate*Gait.Duration); %Time frame of the vedio
Time = ceil(NumberofFrame/StepTime);

Index=1;
%for count=1:StepTime:NumberofFrame-StepTime-DurationTime,    
for Count=StartCount:StepTime:NumberofFrame-DurationTime,      
    RealForeground=read(Gait,Count);
    RealForeground = imresize(RealForeground, 0.5);
    Foreground_BW = rgb2gray(RealForeground);
                  
    %Adjust Intensity Values Using Histogram Equalization
    Foreground_BW = imhistmatch(Foreground_BW, Background_BW);
    Foreground_BW = medfilt2(Foreground_BW,[3 3]);
    %imshow(Foreground_BW);
    ForegroundW= Gait.Width*0.5;
    ForegroundH= Gait.Height*0.5;
   
    %Subtraction Foregroung-Background
    SubtractionF_B1=imabsdiff(Foreground_BW,Background_BW);

    %Apply Median filter to remove Noise
    SubtractionF_B=medfilt2(SubtractionF_B1,[3 3]);
    %imshow(SubtractionF_B);
     
    %Filling black and white
    s=size(SubtractionF_B);
    y=zeros(s(1),s(2));
    for i=1:s(1)
            for j=1:s(2)
                %if(SubtractionF_B(i,j)>Threshold)
                if(SubtractionF_B(i,j)>THRESHOLDSET)
                    FILL_BW(i,j)=0;
                else FILL_BW(i,j)=255;
                end
            end
    end

    [rows, columns]=size(FILL_BW);
    %Convert to Binary Image
    for i=1:rows
        for j=1:columns
            if FILL_BW(i,j) >0
                BinaryImage(i,j)=0;
            else
                BinaryImage(i,j)=1;
            end
        end
    end
    
    %Apply Median filter to remove Noise
    originalImage=medfilt2(BinaryImage,[5 5]);
    
    %Motion Blur
    H = fspecial('motion',MOTIONBLUR,90);
    originalImage = imfilter(originalImage,H,'replicate');

    %Boundary Label the Filtered Image
    [L, num]=bwlabel(originalImage);
    STATS=regionprops(L,'all');
    cc=[];
    removed=0;
    
    %Remove the noisy regions 
    for i=1:num
        dd=STATS(i).Area;
        if (dd < PIXAREA)
            L(L==i)=0;
            removed = removed + 1;
            num=num-1;
        else
        end
    end
    
    [originalImage, num2]=bwlabel(L);
     %imshow(originalImage);
hold on;   

labeledImage = bwlabel(originalImage);

labeledImage = imfill(labeledImage,'holes');

measurements = regionprops(labeledImage, 'BoundingBox', 'Extrema', 'FilledImage', 'Area');


for k = 1 : length(measurements)
  thisBB = measurements(k).BoundingBox;
  rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
  'EdgeColor','g','LineWidth',1 )
end

% Let's extract the second biggest blob - that will be the hand.
allAreas = [measurements.Area];
[sortedAreas, sortingIndexes] = sort(allAreas, 'descend');

try
Body = sortingIndexes(1); % The hand is the second biggest, face is biggest.

catch
 Body = 0;   
  
end

% Use ismember() to extact the hand from the labeled image.
Body = ismember(labeledImage, Body); 

siz=size(labeledImage);
CenterPoint=regionprops(labeledImage,'centroid');
centroids = cat(1, CenterPoint.Centroid);



try
BodyPosition = regionprops(originalImage,'Centroid');
BodyPositionX=BodyPosition.Centroid(:,1); %Body position of X
BodyPositionY=BodyPosition.Centroid(:,2); %Body position of Y

catch

        if (isempty(BodyPosition))
            BodyPositionX=MaxDistance;
            BodyPositionY=PositionX(1);        
        end    
end

%Caculating the distance of the body from the end point.
try
DistanceofBody = BodyPositionX-EndPoint;
catch
    break;
end
PositionofBody=DistanceofBody*OnePixelW;
%fprintf('Position of the long jumper is: x= %0.2f and y= %0.2f\n',BodyPositionX,BodyPositionY);
%fprintf('Distance from the end take-off board is: %0.2f Meters\n',PositionofBody);

BodyP(Index)=PositionofBody;
DistanceD(Index)=RealTotalDistance*(DistanceofBody/DistanceStartToEnd);
RealPositionBodyX(Index)=BodyPositionX;
RealPositionBodyY(Index)=BodyPositionY;
RealFrame(Index)=Count;

% %Create bounding Box
bb=regionprops(labeledImage,'BoundingBox');

% 
% % Crop the individual objects and store them in a cell
n=max(labeledImage(:)); % number of objects
ObjCell=cell(n,1);
for i=1:n
      % Get the bb of the i-th object and offest by 2 pixels in all
      % directions
      bb_i=ceil(bb(i).BoundingBox);
      idx_x=[bb_i(1)-15 bb_i(1)+bb_i(3)+15]; %  bb_i(1) read from 1st location
      idx_y=[bb_i(2)-15 bb_i(2)+bb_i(4)+15];
      if idx_x(1)<1, idx_x(1)=1; end
      if idx_y(1)<1, idx_y(1)=1; end
      if idx_x(2)>siz(2), idx_x(2)=siz(2); end
      if idx_y(2)>siz(1), idx_y(2)=siz(1); end
      % Crop the object and write to ObjCell
      im=labeledImage==i;
      ObjCell{i}=im(idx_y(1):idx_y(2),idx_x(1):idx_x(2));
end

%Plot in the pictures
 P=PositionY(1)*2;
 pos = [EndPoint*2  P;(EndPoint+TenMeters)*2  P;(EndPoint+TenMeters)*2  P;StartPoint*2  P];
 color = {'green','green','green','green'};

 ForegroundM = insertMarker(RealForeground,pos,'x','color',color,'size',8); 
 ForegroundM = insertText(RealForeground, [1 1], sprintf('%0.2f Meters',RealTotalDistance*(DistanceofBody/DistanceStartToEnd)));
 subplot(221);imshow(Background),title('Background');
 subplot(222);imshow(SubtractionF_B1),title('Background Subtraction');
 subplot(223);imshow(ForegroundM),title('Foreground');
 subplot(224);imshow(originalImage),title('Noise Reduction');

%imshow(ForegroundM),title('Foreground'); 

pause(0.01);
  hold on;
  plot(BodyPositionX,BodyPositionY, 'b*');
  %plot(BodyPositionX,BodyPositionY, 'bo');
  %imshow(ObjCell{sortingIndexes(1)});
  hold off;



Index=Index+1;
end
%*********************************************************************

                        %Calculating velocity


%*********************************************************************

for v=1:(Index-DurationTime-1),
    
P1_1=DistanceD(v);

P1_2=DistanceD(v+DurationTime);

rate=(1+(DurationTime*StepTime))*(1/FrameRate);

DistanceA= (P1_1 - P1_2); 

velocityA(v)=DistanceA/rate;

PositionBody(v)=BodyP(v);

fprintf('Velocity at period %d is: %0.2f m/s. at (Distance=%0.3f,Position=%0.3f)\n',v,velocityA(v),DistanceA,PositionBody(v));
end

%**************************************************************************************************************************************************

%Define NaN replace the zero value
inx=1; 
err=0;
Error=1;
si=size(PositionBody);
        %Filtering of Distance position (NaN)
        for inx=1: max(si)-1
         NewPositionBody(inx) = PositionBody(inx);   
         NewvelocityA(inx) = velocityA(inx); 
         NewRealFrame(inx) = RealFrame(inx);
         NewRealPositionBodyX(inx) = RealPositionBodyX(inx);
         NewRealPositionBodyY(inx) = RealPositionBodyY(inx);
         Error(inx)=0;
         
         if(PositionBody(inx)<=0)%-3
         break;
         end
        
        if (velocityA(inx) < 0)% negative speed 
            NewPositionBody(inx) = NaN;
            NewvelocityA(inx) = NaN;
            Error(inx)=1;
        end
         
         if (velocityA(inx) > 13)%maximum speed of a sprinter
            NewPositionBody(inx) = NaN;
            NewvelocityA(inx) = NaN;
            Error(inx)=1;
        end
        
        if (PositionBody(inx)==0 || abs(PositionBody(inx+1)-PositionBody(inx)) >= 5)
            NewPositionBody(inx) = NaN;
            NewvelocityA(inx) = NaN;
            Error(inx)=1;
        end
         
         if (velocityA(inx)==0 || abs(velocityA(inx+1)-velocityA(inx)) >= 5 )
            NewPositionBody(inx) = NaN;
            NewvelocityA(inx) = NaN;
            Error(inx)=1;
         end
        
        inx=inx+1;
        end
        
% Percent of Data Error

ErrorSize=size(Error);
DataError= (sum(Error==1)/max(ErrorSize))*100;
      
%**************************************************************
%Raw Data
%FillPositionBody = PositionBody;
%FillvelocityA = velocityA;

    try
        FillPositionBody = NewPositionBody;
        FillvelocityA = NewvelocityA;
    catch ME
            disp('An error occurred! Please change the background.');
            
            break;
    end

%**************************************************************     
%Fill missing numeric data.
FillPositionBody = Mf(NewPositionBody);    
%Fill missing numeric data.
FillvelocityA = Mf(NewvelocityA); 

x = FillPositionBody; 
y = FillvelocityA; 

FillPositionBody = (x);
FillvelocityA = (y);

%Detection of the startup position
take_inx=1;
start_inx=1;
for take_inx=1:(inx)-1;
    
        if (FillPositionBody(take_inx)<=10)
            break;
        end
    
    start_inx=start_inx+1;
end

%Detection of the end position
take_inx=1;
End_inx=1;
for take_inx=1:(inx)-1;
    
        if (FillPositionBody(take_inx)<=0)
            break;
        end
    
    End_inx=End_inx+1;
end

Startup_pos=FillPositionBody(start_inx);
End_pos=FillPositionBody(End_inx);

Run_distance=Startup_pos;%-End_pos;
Runup_Time=(End_inx*StepTime)/FrameRate;
Runup_Speed=Run_distance/Runup_Time;

%The result of startup position 
fprintf('Data Error(percent): %0.2f\n',DataError);  
fprintf('Time(s): %0.2f \n',Runup_Time);
fprintf('Speed(m/s): %0.2f \n',Runup_Speed);
fprintf('Distance(m): %0.2f \n',Run_distance);

%Ploting velocity graph
plot_inx=1;
for i=start_inx:(End_inx)-1,
 S(plot_inx)= FillPositionBody(i);
 V(plot_inx)= FillvelocityA(i);
 FinalRealFrame(plot_inx) = NewRealFrame(i);
 FinalRealPositionBodyX(plot_inx) = NewRealPositionBodyX(i);
 FinalRealPositionBodyY(plot_inx) = NewRealPositionBodyY(i);
 plot_inx=plot_inx+1;
end

%calculating Vmax

[M,I] = max(smooth(V));
V_Max= V(I);
PositionofV_Max=S(I);
fprintf('Maximum speed is (m/s): %0.2f  \nDistance from the Designated point(m): %0.2f\n' ,V_Max,PositionofV_Max);
VmaxShow=read(Gait,FinalRealFrame(I));
VmaxForeground=read(Gait,FinalRealFrame(I));
VmaxForeground=insertMarker(VmaxForeground,pos,'x','color',color,'size',10); 
VmaxForeground=insertText(VmaxForeground, [FinalRealPositionBodyX(I)*2 FinalRealPositionBodyY(I)*2], sprintf('Vmax:(%0.2f m/s,%0.2f Meters)',V(I),S(I)),'FontSize',18);



%Showing the relationship between Position and Speed
figure1 = figure;
hold all;
% hLine1 = plot ((abs(S-MARKER_CALIBRATION)),(V),'b-*','DisplayName', 'Original data');
hLine1 = plot (smooth(abs(S-MARKER_CALIBRATION)),smooth(V),'b-o','DisplayName', 'Smoothed data');
% hLegend = legend([hLine1,hLine2], 'Location','NorthEast');
hLegend = legend([hLine1], 'Location','NorthEast');
hLegendAxes = axes('Parent',hLegend.Parent, 'Units',hLegend.Units, 'Position',hLegend.Position, ...
                   'XTick',[] ,'YTick',[], 'Color','none', 'YColor','none', 'XColor','none', 'HandleVisibility','off', 'HitTest','off');

hTitle = title(hLegendAxes, 'Data', 'FontWeight','normal', 'FontSize',8);  % Default is bold-11, which is too large

annotation('textbox',...
    [0.15 0.15 0.3 0.10],...
    'String',{['Speed: ',num2str(Runup_Speed) '(m/s) '],['Distance: ',num2str(Run_distance) '(m) '],['Time: ',num2str(Runup_Time), '(s) ']},...
    'FontSize',10,...
    'FontName','Arial',...
    'LineStyle','-',...
    'EdgeColor',[1 1 1],...
    'LineWidth',2,...
    'BackgroundColor',[1.0 1.0 1.0],...
    'Color',[0.0 0.0 0.0]);

xlabel('Position (m)');
ylabel('Speed (m/s)');
title('Relationship between Position and Speed');
%invert position of distance
S=abs(S-MARKER_CALIBRATION);
S=S';

V=V';
%Saving the relationship between Position and Speed
saveas(figure1,(fullfile([datapath,FileName,'MoCap_plot.jpg'])));

%Saving position and speed data
%xlswrite(fullfile([datapath,FileName,'GaitMoCapX.xlsx']),S);
%xlswrite(fullfile([datapath,FileName,'GaitMoCapV.xlsx']),V);
%*********************************************************************************************
%Saving data to the excel file

sheet1 = 'Speed';
%Position and Speed
xlRange = 'A1';
Title_A = {'Speed'};
xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),Title_A,sheet1,xlRange);

xlRange = 'A2';
Label_A = {'Position(m)','Speed (m/s)'};
xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),Label_A,sheet1,xlRange);

xlRange = 'A3';
xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),S,sheet1,xlRange);

xlRange = 'B3';
xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),V,sheet1,xlRange);

%*********************************************************************************************
%Global Distance Time and Speed
xlRange = 'D1';
Title_A = {'Global Distance Time and Speed'};
xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),Title_A,sheet1,xlRange);

xlRange = 'D2';
Label_A = {'Speed (m/s)','Time(S)','Distance(m)'};

xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),Label_A,sheet1,xlRange);

xlRange = 'D3';
xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),Runup_Speed,sheet1,xlRange);

xlRange = 'E3';
xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),Runup_Time,sheet1,xlRange);

xlRange = 'F3';
xlswrite(fullfile([Outputpath,FileName,'Speed.xlsx']),Run_distance,sheet1,xlRange);

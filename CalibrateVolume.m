%Teerawat Kamnardsiri
%Colledge of Arts, Media and Technology, Chiang Mai University
%Gait Speed Detection System
%December, 2022.

clc;
close all;
clear;

datapath = 'S1/';


BGFileName='BG';
BG_File=fullfile([BGFileName,'.MOV']);

%*****************************************************************************************************************
%Read Background Image
BG= VideoReader(sprintf('%s%s',datapath,BG_File));
%Read Background Image
Background = read(BG,2);
Background = imresize(Background, 0.5);
%Size of Background2
Background2W = BG.Width*0.5;
Background2H = BG.Height*.05;

Background_BW = rgb2gray(Background);
Background_BW = imgaussfilt(Background_BW,3);
Background_BW = imadjust(Background_BW,[0.00; 0.4],[0.1; 1.0], 0.25); 
    
%Calibration the system

box_color = {'white'};
Title1 = sprintf(' Callibration the System \n Please perform 2 steps:\n 1: Select area for processing and then double-click.');
Background2 = insertText(Background, [Background2W/2 Background2H/100], Title1, 'FontSize',14,'BoxColor',...
    box_color,'BoxOpacity',0.5,'TextColor','black','AnchorPoint','CenterTop');
figure, imshow(Background2);

h = imrect(gca,[5 100 630 195]); 
position1 = wait(h); % returns coordinates in "position" when user doubleclicks on rectangle

close;

%msgbox('Click OK to continue');  % Wit for user to click OK button.

%Getting the positions of the start point and the end point (3 m)
box_color = {'white'};

Title2 = sprintf(' Steps 2:\n Click mouse at the Marker M1 and then at the Marker M2 (3 m. from M1)\n');
Background2 = insertText(Background, [Background2W/2 Background2H/100], Title2, 'FontSize',14,'BoxColor',...
    box_color,'BoxOpacity',0.5,'TextColor','black','AnchorPoint','CenterTop');
imshow(Background2);
[HorizontalX,HorizontalY]=ginput(2);

close;

Title3 = sprintf(' Callibration for processing have been successfully saved \n');
Background2 = insertText(Background, [Background2W/2 Background2H/100], Title3, 'FontSize',18,'BoxColor',...
    box_color,'BoxOpacity',0.5,'TextColor','black','AnchorPoint','CenterTop');
figure, imshow(Background2);

%Recording calibration positions.
dlmwrite('position1.txt',position1);
type('position1.txt');

dlmwrite('horizontalX.txt',HorizontalX);
type('horizontalX.txt');

dlmwrite('horizontalY.txt',HorizontalY);
type('horizontalY.txt');

fprintf('***************************************************************\n');
fprintf('Callibration data have been successfully saved.');
%pause;

promptMessage = sprintf('Callibration data have been successfully saved.');
button = questdlg(promptMessage, 'OK', 'OK', 'Cancel', 'OK');
if strcmpi(button, 'Cancel')
  return; % Or break or continue
end

close;

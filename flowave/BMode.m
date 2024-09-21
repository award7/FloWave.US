%% BMode2
% Purpose: This code was written to analyze ultrasound (US) B-mode images
% from a digital video recording. The code uses an automated diameter 
% subroutine to determine a blood vessel diameter. 
%
% Inputs: Ultrasound Video (type - '.avi' or '.mov')
% Outputs: Vessel diameter (cm)
%
% Functions: FrameCalibrate, VesselROI, AutoDiameter, imgaussian (written
% by Dirk-Jan Kroon), ginputc (written by Jiro Doke) - See associated
% license agreements for copyright information.
%
% Requires MATLAB R2012 or newer.
%
% Crystal Coolbaugh
% August 4, 2015
% Copyright 2015 Crystal Coolbaugh

%% Format Workspace
clear
format compact
clc
close all
disp('PLEASE FOLLOW THE ON-SCREEN PROMPTS TO PROCESS THE ULTRASOUND DATA.')

%% Identify Folder with Video
% Add video folder to the current path
disp('Choose the file directory that contains the US video.');
DirName=uigetdir;   %Open browser to identify the directory path
addpath(DirName);
ls(DirName); %List file contents for easy copy-past of filename

%% Import US Video Data

% Enter the Video Filename - Can Copy/Paste from List
VideoName = input('Type the video filename and extension (e.g. USVideo.avi): ','s');

%Construct an object using the VideoReader Function
USObj = VideoReader(VideoName);

% ACQUIRE VIDEO SETTINGS: May Vary for Different US or DVD equipment
% Set Image Frame Size 
VidHeight = USObj.Height;
VidWidth = USObj.Width;

% Total Number of Image Frames
nFrames = USObj.NumberofFrames;

% Video Sampling Rate
VidSamRate = round(USObj.FrameRate);

% Video End Time
Duration = USObj.Duration;

%% Pixel Calibration (FrameCalibrate Function)
% Determine the factor needed to convert pixels to the quantity of
% interest: distance. Position of the scales may need to be adjusted for
% different US equipment or screen resolutions. 

%Identify position of the distance scale
OneFrame = read(USObj, 1);
image(OneFrame); colormap gray
title('Define a large ROI around the distance scale (click 1: upper left corner; click 2: lower right corner)');
[X,Y] = ginputc(2, 'Color', 'r', 'LineWidth', 2);
ScaleX = floor(X);
ScaleY = floor(Y);
close all;

ScaleCheck = 0;
disp('Use the cursor to select the scale range (lower & upper extremes).')

while ScaleCheck == 0
    % Set Distance Conversion Factor
    DistScale = OneFrame(ScaleY(1,1):ScaleY(2,1),ScaleX(1,1):ScaleX(2,1));
    DistCon = FrameCalibrate(DistScale);
    
    % Calibration Check
    CalCheck = input('Do you need to repeat the calibration (Y/N)? ','s');
    
    if lower(CalCheck) == 'n'
        ScaleCheck = 1;
    end
end

% set if baseline or FMD measure
mode = input('Baseline or FMD b/f ', 's');

%for FMD enter time offset (when video starts from moment of cuff release)
timeoffset =0;
if mode == 'f'
    timeoffset = input('Enter time offset ' );
end

%set up counters for loop for repeat measures
cont = 'y';
runs =0;

%recalibration rate allows for redesignation of ROI every designated amount
%of seconds
recal = input('Set a recalibration rate? y/n ', 's');
if recal == 'y'
    RedoSec = input('Enter recalibration rate in seconds');
    frameRedo = VidSamRate*RedoSec;
else 
    frameRedo = nFrames+1;
end

%Set times 
    dt = 1/VidSamRate;
    DTime = (0:dt:nFrames*dt-dt)';

    totalDiamsOrig(:,1) = DTime+timeoffset; %will hold raw data
    totalDiamsFilt(:,1) = DTime+timeoffset; %will hold FFT filtered data

%loop to redesignate ROI and recalculate diameter as long as continue = y
while cont == 'y'
    %counter for number of runs
    runs = runs +1;

    %% Calculate the Vessel Diameter
    DiameterPixel = zeros(nFrames,1);

        
    for k = 1:nFrames
        %Read Image at start of Video
        Vessel = read(USObj,k);
        %check if recalibration needs to be done
        redo = mod(k,frameRedo);
        if k==1 || redo == 0
            %Select ROI around the vessel 
            newROI ='y';
            while newROI == 'y'
            figure; imagesc(Vessel);
            [DROIx,DROIy,Center,theta] = VesselROI(Vessel);
             % Restrict image area to ROI
            DiameterImage = Vessel(DROIy(1,1):DROIy(2,1),...
                DROIx(1,1):DROIx(2,1));
            AutoDiameter(DiameterImage,Center',theta(1),1);
            
            newROI = input('Change ROI? y/n ', 's');

            end
        end

      
        % Restrict image area to ROI
        DiameterImage = Vessel(DROIy(1,1):DROIy(2,1),...
                DROIx(1,1):DROIx(2,1));
        %Calculate the Diameter
        DiameterPixel(k,1) = AutoDiameter(DiameterImage,Center',theta(1),0);
   
    end
    
        


    % Replace Gaps and Errant Values with Mean Diameter
    MeanD = mean(DiameterPixel(intersect(find(DiameterPixel~= 200),find(DiameterPixel~=0))));
    DiffD = diff(DiameterPixel(find(DiameterPixel~=0)))/MeanD;

    %counter for frames that did not return a value
    nReplacedNaN = 0;
    for i = 1:length(DiameterPixel)
        if DiameterPixel(i,1) == 200
            DiameterPixel(i,1) = MeanD;
            nReplacedNaN = nReplacedNaN+1;
        end
    end
    
    % If Diameter Varies by > 30% from the Mean, Replace with Mean Value
    % Counter for replaced outliers
    nReplacedOut = 0;
    for i = 1:length(DiffD)
        if abs(DiffD(i,1)) > .3
            DiameterPixel(i,1) = MeanD;
            nReplacedOut = nReplacedOut+1;
        end
    end
    
    %percent of frames read
    percentReadFrames = (nFrames- nReplacedNaN)/nFrames*100

    % Convert Diameter to Centimeters
    Diameter = DiameterPixel.*DistCon;
    MeanDiam = mean(Diameter);
    StdDiam = std(Diameter);


%FMD noise reduction and data smoothing 
if mode =='f'
   totalDiamsOrig(:, runs+1) = Diameter;
   
   %FastFourierTransformFilter
   [DFilt, a,b,c]= fftf(DTime, Diameter, [],150);
   DFilt=transpose(DFilt);
   totalDiamsFilt(:,runs+1)= DFilt;
    
   %averageFFT filter over 2 second intervals
   pointsFFT = VidSamRate*2;
   for i = 1:((length(Diameter)/(pointsFFT)))-1
       totalDiamsCombo(i,1) = mean(DTime((((i-1)*pointsFFT)+1):i*pointsFFT,:))+timeoffset;
       totalDiamsCombo(i,runs+1) = mean(DFilt((((i-1)*pointsFFT)+1):i*pointsFFT,:));
   end
   
   %Find maxes of FFT & 2 sec average
   [maxDCombo, maxTComboPos] = max(totalDiamsCombo(:,runs+1));
   maxTCombo = totalDiamsCombo(maxTComboPos,1);
   maxDiam(1,runs)=maxDCombo;
   maxDiam(2,runs)=maxTCombo;
   maxDiam

else
     % Filter Diameter - Smooth with Savitsky Golay Filter
    MeanDiam
    DFilt = sgolayfilt(Diameter,3,11);
    totalDiamsOrig(:, runs+1) = Diameter;
     
end
    
%% Plot and Summarize

% Error Limits - +/- 10% from Mean Diameter
% User Inspection of Data - If diameter exceeds error limits, repeat
% epoch selection or enter a single diameter value.
for i = 1:length(totalDiamsOrig)
    PosError(i,1) = mean(totalDiamsOrig(:,runs+1)) + .1*mean(totalDiamsOrig(:,runs+1));
    NegError(i,1) = mean(totalDiamsOrig(:,runs+1)) - .1*mean(totalDiamsOrig(:,runs+1));
end

%original
figure, grid on, hold on;
plot(DTime+timeoffset, totalDiamsOrig(:,runs+1),'b','DisplayName','Original Diameter', 'LineWidth', 2); hold on;
plot(DTime+timeoffset,PosError,'k','LineWidth',2); hold on;
plot(DTime+timeoffset,NegError,'k', 'LineWidth',2); hold on;
title('Filtered Vessel Diameter & +/- 10% Error Limits');
xlabel('Time(s)');
ylabel('Diameter (cm)');
legend('Diameter');

if mode == 'f'
    plot(totalDiamsCombo(:,1),totalDiamsCombo(:,runs+1),'r','DisplayName','Smoothed Diameter', 'LineWidth', 2); hold on;
    legend('Raw Diameter', 'Filtered Diameter');
    hold off;
end

cont = input('Continue? y/n ', 's');
end



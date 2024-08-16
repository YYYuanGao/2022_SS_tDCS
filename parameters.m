% BY Ke Jia: Version 1_20220623 15:31
% list all the parameters used in this experiment

Param = struct;

%% Screen Settings
Param.Settings.ViewDistance      = 720;              % 1100 mm 
Param.Settings.ScrnResolution    = [0 0 1920 1080];   % rect_computer = [0 0 40 30];
Param.Settings.SquareSize        = 1920;              % 400 mm
Param.Settings.SquareLength      = 532;               % 154 mm 
Param.Settings.PixelPerDegree    = 2*Param.Settings.ViewDistance*tan(1/2/180*pi)*Param.Settings.SquareSize/Param.Settings.SquareLength;       

%% Keys for response
Param.Keys.Space    = 32;  
Param.Keys.Esc      = 27;
Param.Keys.Left     = 37;  
Param.Keys.Right    = 39;
Param.Keys.Down     = 40;
Param.Keys.Trigger1 = 83;  % 's'       

%% parameters for stimuli
% Locations
Param.Stimuli.Eccentricity   = 5; % degree
Param.Stimuli.Locations(1,:) = [Param.Settings.ScrnResolution(3)/2 - Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(2,:) = [Param.Settings.ScrnResolution(3)/2 + Param.Stimuli.Eccentricity * Param.Settings.PixelPerDegree,Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.Locations(3,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2] + sqrt(2)/2*Param.Stimuli.Eccentricity*Param.Settings.PixelPerDegree;
Param.Stimuli.Locations(4,:) = [Param.Settings.ScrnResolution(3)/2 , Param.Settings.ScrnResolution(4)/2];
Param.Stimuli.LocationsText  = {'Left','Right','Lower Right','Centre'};

% Gratings
Param.Stimuli.OuterRadius      = 2.5;
% Param.Stimuli.OuterRadius_edge = Param.Stimuli.OuterRadius - 0.5;
Param.Stimuli.InnerRadius      = 1;
Param.Stimuli.OuterSize        = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
Param.Stimuli.InnerSize        = round(Param.Stimuli.InnerRadius*Param.Settings.PixelPerDegree);   % radius pixels
Param.Stimuli.Inneroffset      = 2; 

Param.Stimuli.GratingOri       = 45; 
Param.Stimuli.ReferenceOri     = 60; 
Param.Stimuli.OriJitter        = 5; 

Param.Stimuli.GratingContrast  = 0.7;
Param.Stimuli.Spatial_freq     = 2.5/Param.Settings.PixelPerDegree;       % 0.02 % 6 bars in total
Param.Stimuli.SmoothSD         = Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree/3;          % 0.5degree*1/3
 
%% parameters for fixation
Param.Fixation.CrossColor    = [1,1,1]*255;
Param.Fixation.CrossSize     = 0.3*Param.Settings.PixelPerDegree;
Param.Fixation.CrossWidth    = 0.1*Param.Settings.PixelPerDegree;
% +
Param.Fixation.CrossLoc      = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize, Param.Settings.ScrnResolution(3)/2, Param.Settings.ScrnResolution(3)/2;...
                                     Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize] + [offset',offset',offset',offset'];
% X
Param.Fixation.CrossLoc2     = [Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(3)/2+Param.Fixation.CrossSize*sqrt(2)/2;...
                                     Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2+Param.Fixation.CrossSize*sqrt(2)/2, Param.Settings.ScrnResolution(4)/2-Param.Fixation.CrossSize*sqrt(2)/2]+[offset',offset',offset',offset'];  
Param.Fixation.OvalSize      = 2*Param.Settings.PixelPerDegree;
Param.Fixation.OvalColor     = [0,0,0]; 

%% parameters for trials
Param.Trial.ITI              = 0.5;     
Param.Trial.ISI              = 0.5;    % 12 frames= 600 ms
Param.Trial.StiDura          = 0.4; 
Param.Trial.MaxRT            = 1.5; 
Param.Trial.Duration         = 3;

%% parameters for prac-test 
Param.Prac.TrialNum = 30;
Param.Prac.DeltaAngle = 10;

%% parameters for testing_method of constant stimuli 
Param.conTest.TrialNum_perdiff = 20;
Param.conTest.Taskdiff         = [3 2 1 0 -1 -2 -3]*5;
Param.conTest.TrialNum         = Param.conTest.TrialNum_perdiff * length(Param.conTest.Taskdiff); 
Param.conTest.TrialNum_minirun = 50;
%% parameters for testing_staircase
Param.stairTest.start = 10;
Param.stairTest.MaxTrial = 60;
Param.stairTest.Up = 1;     %increase after 1 wrong
Param.stairTest.Down = 2;   %decrease after 3 consecutive right
Param.stairTest.StepSizeDown = 0.5;         
Param.stairTest.StepSizeUp = 0.5;
Param.stairTest.StopCriterion = 'reversals';   
Param.stairTest.StopRule = 10;
Param.stairTest.ReversalsUnused = 4; 
Param.stairTest.xMax = 20;
Param.stairTest.xMin = 0;

%% open window 
screens = Screen('Screens');
screenNumber = max(screens);	
Screen('Preference', 'SkipSyncTests', 0);

white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
gray = round((white+black)/2);
wnd = Screen('OpenWindow',screenNumber,gray,Param.Settings.ScrnResolution); 
Screen('BlendFunction',wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

RefreshDur = Screen('GetFlipInterval',wnd);
Slack = RefreshDur/2;

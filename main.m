% SS_main: this calls the main functions of the tDCS surround suppresion project
% subjID: e.g. S001

clear all;
clc;

SubjID  = '111';
sess    = 1;
run     = 1;
run_case= 2;

%%
offset = [0,0];
CurrDir = pwd;
SetupRand;
% set_test_gamma;
HideCursor;
parameters;

%%

switch run_case 
    case 0
        sample;  
      
    case 1
       % prac       
        location_used = 3; 
        Param.Stimuli.Spatial_freq  = 1/Param.Settings.PixelPerDegree;
        Param.Stimuli.OuterRadius   = 1.5;
        Param.Stimuli.OuterSize     = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
        ori_prac; 
     
    case 2   
        % stair       
        location_used = 3; 
        Param.Stimuli.Spatial_freq  = 1/Param.Settings.PixelPerDegree;
        Param.Stimuli.OuterRadius   = 1.5;
        Param.Stimuli.OuterSize     = round(Param.Stimuli.OuterRadius*Param.Settings.PixelPerDegree);   % radius pixels
        ori_test_staircase;  
    
    case 3      
        location_used = 4; 
        Param.Stimuli.Spatial_freq  = 2.25/Param.Settings.PixelPerDegree;
        SS_constant;
end

delete *.asv

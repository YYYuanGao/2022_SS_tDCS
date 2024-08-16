
% Inner sti
[x,y]=meshgrid(-Param.Stimuli.InnerSize:Param.Stimuli.InnerSize,-Param.Stimuli.InnerSize:Param.Stimuli.InnerSize);
phase = rand*2*pi; 
angle = 45/180*pi;
sti1_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
sti1_final(sqrt(x.^2 + y.^2) > Param.Stimuli.InnerSize-Param.Stimuli.Inneroffset) = gray; %circle mask

mm = Screen('MakeTexture', wnd, sti1_final); 
Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(3,1)-Param.Stimuli.InnerSize,Param.Stimuli.Locations(3,2)-Param.Stimuli.InnerSize,Param.Stimuli.Locations(3,1)+Param.Stimuli.InnerSize,Param.Stimuli.Locations(3,2)+Param.Stimuli.InnerSize]);

% Outer sti
[x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
phase = rand*2*pi; 
angle = 60/180*pi;
sti2_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
sti2_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
sti2_final(sqrt(x.^2 + y.^2) < Param.Stimuli.InnerSize) = gray; %circle mask

annulus=ones(Param.Stimuli.OuterSize*2+1)*white;
annulus(sqrt(x.^2 + y.^2) < Param.Stimuli.InnerSize) = black; %circle mask 
sti_final = repmat(sti2_final,[1,1,3]); 
sti_final(:,:,4) = annulus;

mm = Screen('MakeTexture', wnd, sti_final);   
Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(3,1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(3,2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(3,1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(3,2)+Param.Stimuli.OuterSize]);

% fixation
Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
vbl = Screen('Flip',wnd);

     
%% close all
is_true = 0;
while ~is_true
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.Esc)
        is_true = 1;
    end
end

Screen('CloseAll');
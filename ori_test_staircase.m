% orientation discrimination task
% staircase 2down1up
%  By Ke Jia 2022/06/24 09:02
% modified by Yuan Gao 

%% check parameters used
if exist([CurrDir '\Results\staircaseTest\' SubjID '\' SubjID '_test_results_sess' num2str(sess) '_run' num2str(run) '.mat'],'file')
    disp(' ');
    disp([SubjID '_run' num2str(run) ' has been test, please enter a new run number!']); 
    disp(' ');
    abort;      
end

results = zeros(Param.stairTest.MaxTrial,15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results:
% 1 trial number              %10    trial onset
% 2 sti_loc                   %11    fix onset
% 3 sti_ori                   %12    ITI
% 4 SN of st1                 %13    st1
% 5 SN of st2                 %14    ISI
% 6 response                  %15    st2
% 7 task_diff
% 8 accuracy
% 9 reaction time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% staircase settings
task_dif = Param.stairTest.start;
UD = PAL_AMUD_setupUD('up',Param.stairTest.Up,'down',Param.stairTest.Down);
UD = PAL_AMUD_setupUD(UD,'StepSizeDown',Param.stairTest.StepSizeDown,'StepSizeUp',Param.stairTest.StepSizeUp,'stopcriterion',Param.stairTest.StopCriterion,'stoprule',Param.stairTest.StopRule,'startvalue',task_dif,...
    'xMax',Param.stairTest.xMax,'xMin',Param.stairTest.xMin,'truncate','yes');

%% start
DrawFormattedText(wnd,'Press space to start!','center','center', black);
Screen('Flip',wnd);
is_true = 0;
while (is_true == 0)
    [~,~,keyCode] = KbCheck;
    if keyCode(Param.Keys.Space)
        is_true = 1;
    elseif keyCode(Param.Keys.Esc)
        abort;
    end
end

%% Go!
while (~UD.stop & (size(UD.response)<Param.stairTest.MaxTrial))
    trial_i = size(UD.response,2)+1;
    results(trial_i,1) = trial_i; 
    results(trial_i,7) = UD.xCurrent;
       
    results(trial_i,2) = location_used;
    results(trial_i,3) = Param.Stimuli.GratingOri;
      
    trial_onset = GetSecs;
    results(trial_i,10) = trial_onset; 

    trial_delta_ori= sign(rand-0.5); 
    if trial_delta_ori == 0
        trial_delta_ori = 1;
    end

    trial_order = ((rand-0.5)>=0);          
    if trial_order    %1ref first %0 test first
        results(trial_i,4) = results(trial_i,3) + (rand-0.5)*2*Param.Stimuli.OriJitter;
        results(trial_i,5) = results(trial_i,4) + trial_delta_ori*results(trial_i,7);
    else
        results(trial_i,5) = results(trial_i,3) + (rand-0.5)*2*Param.Stimuli.OriJitter;
        results(trial_i,4) = results(trial_i,5) + trial_delta_ori*results(trial_i,7);    
    end

    %% ITI
    Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd);
    results(trial_i,11) = vbl-trial_onset; 
  
    %% sti1
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = results(trial_i,4)/180*pi;
    sti1_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti1_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
    
    mm = Screen('MakeTexture', wnd, sti1_final); 
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI); %-Slack
    results(trial_i,12) = vbl-trial_onset-results(trial_i,11); 

    %% ISI
    Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura); %-Slack
    results(trial_i,13) = vbl-trial_onset-sum(results(trial_i,11:12)); 

    %% sti2
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = results(trial_i,5)/180*pi;
    sti2_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti2_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
    
    mm = Screen('MakeTexture', wnd, sti2_final); 
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI); %-Slack
    stimulus_onset = vbl;
    results(trial_i,14) = vbl-trial_onset-sum(results(trial_i,11:13));

    %% response collection
    Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);%-Slack
    results(trial_i,15) = vbl-trial_onset-sum(results(trial_i,11:14));      
    
    is_true = 0;
    while (is_true == 0 & GetSecs - stimulus_onset < Param.Trial.MaxRT)
        [~,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right)
            results(trial_i,6) = 2;              %6 clockwise
            if results(trial_i,5)> results(trial_i,4)
                results(trial_i,8) = 1;          %8 correct or not
            end
            results(trial_i,9) = RT_time - stimulus_onset;  %9 RT
            is_true = 1;
        elseif keyCode(Param.Keys.Left)
            results(trial_i,6) = 1;              
            if results(trial_i,5) < results(trial_i,4)
                results(trial_i,8) = 1;
            end
            results(trial_i,9) = RT_time - stimulus_onset;
            is_true = 1;
        elseif keyCode(Param.Keys.Esc)
            abort;
        end
    end

    UD = PAL_AMUD_updateUD(UD, results(trial_i,8)); %update UD structure
           
    while (GetSecs - trial_onset < Param.Trial.Duration)
    end         
end
       
%% save the data
threshold_value = PAL_AMUD_analyzeUD(UD, 'reversals', max(UD.reversal)-Param.stairTest.ReversalsUnused);
resultsDir = [CurrDir '\Results\staircaseTest\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder和mkder可以用来创建新文件夹
    mkdir(resultsDir);
end

cd(resultsDir);
results_name = [SubjID '_test_results_sess' num2str(sess) '_run' num2str(run) '.mat'];
save(results_name,'results','UD','threshold_value','Param');
cd(CurrDir); 

%% plot
end_trial = find(UD.reversal==max(UD.reversal));  
task_diff_temp = abs(results(:,4) - results(:,5));
figure;
plot(1:end_trial,task_diff_temp(1:end_trial));
axis([0 Param.stairTest.MaxTrial 0 10]);  

%%
% reset_test_gamma;
ShowCursor;
Screen('CloseAll');
 
delete *.asv

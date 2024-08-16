% surround supression task---method of constant stimuli
% By Ke Jia 2022/06/24 09:36
% modified by Yuan Gao

%% check parameters used
if exist([CurrDir '\Results\constantTest\' SubjID '\' SubjID '_test_results_sess' num2str(sess) '_run' num2str(run) '.mat'],'file')
    disp(' ');
    disp([SubjID '_run' num2str(run) ' has been test, please enter a new run number!']); 
    disp(' ');
    abort;      
end

results = zeros(Param.conTest.TrialNum,17);
trial_index = randperm(Param.conTest.TrialNum);
trial_index = mod(trial_index,length(Param.conTest.Taskdiff));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results:
% 1 trial number              %10    trial onset
% 2 sti_loc                   %11    fix onset
% 3 sti_ori                   %12    ITI
% 4 ori of st1                %13    st1
% 5 ori of st2                %14    ISI
% 6 response                  %15    st2
% 7 task_diff                 %16    ori of large grating
% 8 accuracy                  %17    response relative to ref: if clockwise
% 9 reaction time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Go!
for trial_i = 1:Param.conTest.TrialNum

    % start
    if mod(trial_i,Param.conTest.TrialNum_minirun)==1  
        
        if trial_i ~= 1  
            DrawFormattedText(wnd,'Take a rest! Press space to start!','center','center', black);
        else
            DrawFormattedText(wnd,'Press space to start!','center','center', black);
        end
        
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
    end
    
    results(trial_i,1) = trial_i; 
    
    Curr_cond = trial_index(trial_i);
    if Curr_cond == 0 
        Curr_cond = length(Param.conTest.Taskdiff); 
    end
   
    results(trial_i,7) = Param.conTest.Taskdiff(Curr_cond);  
    results(trial_i,2) = location_used;
    results(trial_i,3) = Param.Stimuli.GratingOri;
      
    trial_onset = GetSecs;
    results(trial_i,10) = trial_onset; 

    trial_order = ((rand-0.5)>=0);  
    ori_temp = (rand-0.5)*2*Param.Stimuli.OriJitter;
    results(trial_i,16)= Param.Stimuli.ReferenceOri + ori_temp;
    if trial_order    %1ref first %0 test first        
        results(trial_i,4) = results(trial_i,3) + ori_temp;
        results(trial_i,5) = results(trial_i,4) + results(trial_i,7);
    else
        results(trial_i,5) = results(trial_i,3) + ori_temp;
        results(trial_i,4) = results(trial_i,5) + results(trial_i,7);
    end


    %% ITI
    Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd);
    results(trial_i,11) = vbl-trial_onset; 
    
    
    %% sti1
    [x,y]=meshgrid(-Param.Stimuli.InnerSize:Param.Stimuli.InnerSize,-Param.Stimuli.InnerSize:Param.Stimuli.InnerSize);
    phase = rand*2*pi; 
    angle = results(trial_i,4)/180*pi;
    sti1_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti1_final(sqrt(x.^2 + y.^2) > Param.Stimuli.InnerSize-Param.Stimuli.Inneroffset) = gray; %circle mask
    
    mm = Screen('MakeTexture', wnd, sti1_final); 
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.InnerSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.InnerSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.InnerSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.InnerSize]);
%     Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
%     Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
       
    if trial_order
        [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
        phase = rand*2*pi; 
        angle = results(trial_i,16)/180*pi;
        sti3_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
        sti3_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
        
        annulus=ones(Param.Stimuli.OuterSize*2+1)*white;
        annulus(sqrt(x.^2 + y.^2) < Param.Stimuli.InnerSize) = black; %circle mask
        sti_final = repmat(sti3_final,[1,1,3]); 
        sti_final(:,:,4) = annulus;

        mm = Screen('MakeTexture', wnd, sti_final);   
        Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    end

    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI); %-Slack
    results(trial_i,12) = vbl-trial_onset-results(trial_i,11); 

    %% ISI
%     Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
%     Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    DrawFormattedText(wnd,'','center','center', black);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura); %-Slack
    results(trial_i,13) = vbl-trial_onset-sum(results(trial_i,11:12)); 

    %% sti2
    [x,y]=meshgrid(-Param.Stimuli.InnerSize:Param.Stimuli.InnerSize,-Param.Stimuli.InnerSize:Param.Stimuli.InnerSize);
    phase = rand*2*pi; 
    angle = results(trial_i,5)/180*pi;
    sti2_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti2_final(sqrt(x.^2 + y.^2) > Param.Stimuli.InnerSize-Param.Stimuli.Inneroffset) = gray; %circle mask
    
    mm = Screen('MakeTexture', wnd, sti2_final); 
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.InnerSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.InnerSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.InnerSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.InnerSize]);
%     Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
%     Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);
        
    if ~trial_order
        [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
        phase = rand*2*pi; 
        angle = results(trial_i,16)/180*pi;
        sti3_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
        sti3_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
        
        annulus=ones(Param.Stimuli.OuterSize*2+1)*white;
        annulus(sqrt(x.^2 + y.^2) < Param.Stimuli.InnerSize) = black; %circle mask
        sti_final = repmat(sti3_final,[1,1,3]); 
        sti_final(:,:,4) = annulus;

        mm = Screen('MakeTexture', wnd, sti_final);   
        Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    end

    vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI); %-Slack
    stimulus_onset = vbl;
    results(trial_i,14) = vbl-trial_onset-sum(results(trial_i,11:13));

    %% response collection
%     Screen(wnd,'FillOval', gray, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
%     Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    DrawFormattedText(wnd,'','center','center', black);
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);     %-Slack   
    results(trial_i,15) = vbl-trial_onset-sum(results(trial_i,11:14));

    is_true = 0;
    while (is_true == 0 & GetSecs - stimulus_onset < Param.Trial.MaxRT)
        [~,RT_time,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right)
            results(trial_i,6) = 2;              %6 sti in interval 2
            if results(trial_i,5)> results(trial_i,4)
                results(trial_i,8) = 1;          %8 被试反应正确与否
            end
            results(trial_i,9) = RT_time - stimulus_onset;  %9 反应时
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

    if trial_order %1ref first %0 test first 
        results(trial_i,17) =  results(trial_i,6)-1;
    else
        results(trial_i,17) =  2-results(trial_i,6);
    end

           
    while (GetSecs - trial_onset < Param.Trial.Duration)
    end         
end
       
%% fit the data
StimLevels = Param.conTest.Taskdiff;
OutOfNum = ones(1,length(Param.conTest.Taskdiff))*Param.conTest.TrialNum_perdiff;
PF = @PAL_Logistic;
paramsValues = [0 1 0 0];
paramsFree = [1 1 0 0];

figure;
% calculate the accuracy
NumPos = zeros(1,length(Param.conTest.Taskdiff));
for diff_i = 1:length(Param.conTest.Taskdiff)
    NumPos(1,diff_i) = sum(results((results(:,7) == Param.conTest.Taskdiff(diff_i)),17));
end

[paramsValues,a,b,c] = PAL_PFML_Fit(StimLevels,NumPos,OutOfNum,paramsValues,paramsFree,PF);

PropCorrectData = NumPos./OutOfNum;  
StimLevelsFine = [min(StimLevels):(max(StimLevels)-min(StimLevels))./1000:max(StimLevels)];
Fit = PF(paramsValues,StimLevelsFine);
plot(StimLevels,PropCorrectData,'k.' ,'markersize' ,40);
set(gca,'fontsize',12);
axis([-Param.conTest.Taskdiff(1) Param.conTest.Taskdiff(1) 0 1]);
hold on;
plot(StimLevelsFine,Fit,'g-' ,'linewidth' ,4);

%% save the data
resultsDir = [CurrDir '\Results\constantTest\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder和mkder可以用来创建新文件夹
    mkdir(resultsDir);
end

cd(resultsDir);
results_name = [SubjID '_test_results_sess' num2str(sess) '_run' num2str(run) '.mat'];
save(results_name,'results','paramsValues','Param');
cd(CurrDir); 

%%
% reset_test_gamma;
ShowCursor;
Screen('CloseAll');
 
delete *.asv

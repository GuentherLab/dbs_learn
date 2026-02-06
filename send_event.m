function [event_times,code_sent,dio] = send_event(eventCode,dio,Interval,Dur,Rep,devID,varargin)
    %% output
    % 1) "event_times": [N x 1]
    %% event code
    % # "eventCode" =  chans to send envent, e.g. [1,4] => simultaneous DO-1 and DO-4
    % # "eventCodeVect" = converts "eventCode" to 1x8 binary vector; e.g. [1,4] => [1,0,0,1,0,0,0,0]
    eventCodeVect = false(1,8);
    if isempty(eventCode)
        eventCode = 0;
    else
        % convert coode
        assert(all((eventCode<=8) & (eventCode>=1)),'"eventCode" must be a vector of integer elements from 1 to 8')
        eventCodeVect(eventCode) = true;
    end
    
    %%
    default_Interval = 0.5; % 0.5 sec
    default_Duration = 0.05; % 0.1 sec
    default_Repeat = 1; % 4
    default_devID = 'Dev1'; % USB6009
    default_verbose = 1;
    %
    if nargin<2, dio = []; end  
    if nargin<3, Interval = default_Interval; end % interval = dur + delay
    if nargin<4, Dur = default_Duration; end
    if nargin<5, Rep = default_Repeat; end
    if nargin<6, devID = default_devID; end
    if nargin>=7, verbose = varargin{1}; else; verbose = default_verbose; end
    assert(Dur<Interval,'"Dur" must be less than "Interval"')
    %}
    %%
    
    % establish dio session
    if isempty(dio) || ~isvalid(dio)
        %global dio % this is not idea but let's keep it this way
        dio = digitalio('nidaq',devID);
        addline(dio, 0:7, 'out');
    end
    event_times = [];
    t = timer('StartFcn',@(~,~)fprintf('timer started, %1.0f beacons in 0.5 sec.\n',Rep),...
          'StartDelay',0.5,...
          'TimerFcn',@(~,~) toggle,...
          'Period',Interval,...
          'TasksToExecute',Rep,...
          'ExecutionMode','fixedRate',...
          'StopFcn',@TimerCleanup);
    start(t)
    wait(t)
    delete(t)
    code_sent = repmat(eventCodeVect,Rep,1);
    
    function toggle(~,~)
            if verbose, fprintf('Beacon On\n'); end
            % eventCodeVect = [0 1 0 0 0 0 0 0]; % beacon only
            %eventCodeVect = [1 1 0 0 0 0 1 1]; % all 
            %eventCodeVect = [1 1 0 0 0 0 1 1]; % screenclicker + iDBS 
            %         eventCodeVect = [1 1 1 1 1 1 1 1]
            event_times = cat(1,event_times,now());
            putvalue(dio, eventCodeVect);
            pause(Dur)
            resetCodeVect = false(1,8);
            putvalue(dio,resetCodeVect);
            if verbose, fprintf('Beacon Off\n'); end

    end

    function TimerCleanup(mTimer,~)
        if verbose, disp('Stopping Timer.'); end
        %delete(mTimer)
    end
    % # notes:
    %   1) ACTIVA @ a0.6, a0.7 (idx 7,8)
    %   2) screenClicker @ a.0.0 (idx 1)
    %   3) beacon @ a.0.1 (idx 2)
    
    %clear global dio

end
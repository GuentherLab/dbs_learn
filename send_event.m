function [event_times,code_sent,dio] = send_event(eventCode,varargin)
    %% output
    % 1) "event_times": [N x 1]
    %% event code
    % # "eventCode" =  chans to send envent, e.g. [1,4] => simultaneous DO-1 and DO-4
    % # "eventCodeVect" = converts "eventCode" to 1x8 binary vector; e.g. [1,4] => [1,0,0,1,0,0,0,0]
    tic1 = tic;
    eventCodeVect = false(1,8);
    if isempty(eventCode)
        eventCode = 0;
    else
        % convert coode
        assert(all((eventCode<=8) & (eventCode>=1)),'"eventCode" must be a vector of integer elements from 1 to 8')
        eventCodeVect(eventCode) = true;
    end
    
    %%
    p = inputParser; % still faster
    addParameter(p, 'dio',         []);      % existing session (or [] to create)
    addParameter(p, 'Interval',    0.5);    % interval between beacons (sec)
    addParameter(p, 'Dur',         0.05);   % pulse duration (sec)
    addParameter(p, 'Rep',         1);      % number of repeats
    addParameter(p, 'DevID',       'Dev1'); % NI device ID (USB6009)
    addParameter(p, 'StartDelay', 0);      % delay before first beacon (sec)
    addParameter(p, 'verbose',     1);
    addParameter(p, 'debug_timer', false); % bypass dio; test timer timing only
    parse(p, varargin{:});
    dio         = p.Results.dio;
    Interval    = p.Results.Interval;
    Dur         = p.Results.Dur;
    Rep         = p.Results.Rep;
    devID       = p.Results.DevID;
    start_delay = p.Results.StartDelay;
    verbose     = p.Results.verbose;
    debug_timer = p.Results.debug_timer;
    if Rep > 1
        assert(Dur<Interval,'"Dur" must be less than "Interval"')
    end
    %}
    
    %%
    
    % establish dio session
    if ~debug_timer
        if isempty(dio) || ~isvalid(dio)
            %global dio % this is not idea but let's keep it this way
            dio = digitalio('nidaq',devID);
            addline(dio, 0:7, 'out');
            data = [0 0 0 0 0 0 0 0];
            putvalue(dio,data);
        end
    else
        if verbose, fprintf('[debug_timer] dio connection skipped.\n'); end
    end
    event_times = [];
    if Rep == 1
        % fast path: skip timer object overhead for single pulse
        if verbose, fprintf('single pulse, StartDelay=%1.3f sec.\n', start_delay); end
        if start_delay > 0, pause(start_delay); end
        
        toggle();
    else
        t = timer('StartFcn',@(~,~)fprintf('timer started, %1.0f beacons in %1.3f sec.\n',Rep,start_delay),...
              'StartDelay',start_delay,...
              'TimerFcn',@(~,~) toggle,...
              'Period',Interval,...
              'TasksToExecute',Rep,...
              'ExecutionMode','fixedRate',...
              'StopFcn',@TimerCleanup);
        start(t)
        if isvalid(t), wait(t); end
        if isvalid(t), delete(t); end
    end
    toc(tic1)
    code_sent = repmat(eventCodeVect,Rep,1);
    
    function toggle(~,~)
            if verbose, fprintf('Beacon On\n'); end
            % eventCodeVect = [0 1 0 0 0 0 0 0]; % beacon only 
            %eventCodeVect = [1 1 0 0 0 0 1 1]; % all 
            %eventCodeVect = [1 1 0 0 0 0 1 1]; % screenclicker + iDBS 
            %         eventCodeVect = [1 1 1 1 1 1 1 1]
            event_times = cat(1,event_times,now());
            if ~debug_timer
                putvalue(dio, eventCodeVect);
                pause(Dur)
                resetCodeVect = false(1,8);
                putvalue(dio,resetCodeVect);
            else
                pause(Dur)
                fprintf('[debug_timer] putvalue skipped (Dur=%.3fs).\n', Dur);
            end
            if verbose, fprintf('Beacon Off\n'); end

    end

    function TimerCleanup(mTimer,~)
        if verbose, disp('Stopping Timer.'); end
        delete(mTimer)
    end
    % # notes:
    %   1) ACTIVA @ a0.6, a0.7 (idx 7,8)
    %   2) screenClicker @ a.0.0 (idx 1)
    %   3) beacon @ a.0.1 (idx 2)
    
    %clear global dio

end
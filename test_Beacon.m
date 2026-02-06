function beacon_times = test_Beacon(varargin)
    Interval =0.5; % 0.5 sec
    Dur = 0.05; % 0.1 sec
    Rep = 4; % 4
    
    if nargin>=1, Interval = varargin{1}; end % interval = dur + delay
    if nargin>=2, Dur = varargin{2}; end
    if nargin>=3, Rep = varargin{3}; end
    if nargin>=4, verbose = varargin{4}; else; verbose = 1; end
    assert(Dur<Interval,'"Dur" must be less than "Interval"')
    %%
    global dio % this is not idea but let's keep it this way
    
    % establish dio session
    if isempty(dio) || ~isvalid(dio)
        dio = digitalio('nidaq','Dev1');
        addline(dio, 0:7, 'out');
    end
    beacon_times = [];
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
    
    
    function toggle(~,~)
            if verbose, fprintf('Beacon On\n'); end
            beacon_times = cat(2,beacon_times,now());

            data = [0 1 0 0 0 0 0 0]; % beacon only
            %data = [1 1 0 0 0 0 1 1]; % all 
            %data = [1 1 0 0 0 0 1 1]; % screenclicker + iDBS 
            
            %         data = [1 1 1 1 1 1 1 1]
            putvalue(dio, data);
            pause(Dur)
            data = [0 0 0 0 0 0 0 0];
            putvalue(dio,data);
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
    
    clear global dio

end
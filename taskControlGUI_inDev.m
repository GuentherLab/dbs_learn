function figTC=taskControlGUI(taskState,pulseParam,beacon_times_filepath,beacon_times)

    pause_requested = taskState.pause_requested;
    pause_isActive = taskState.pause_isActive;
    task_isRunning = taskState.task_isRunning;

    % Initialize global variables if they do not exist
    if isempty(task_isRunning) || ~task_isRunning
        warning('No active FLvoice task running...')
        figTC = [];
        return
    end
    stop_requested = false;
    if ~isfield(taskState,'stop_requested'), taskState.stop_requested = false; end

    % Create a figure for the GUI
    figTC = figure('Name', 'Task Control', 'NumberTitle', 'off', ...
                 'Position', [100, 100, 420, 200], 'Tag', 'TaskControlGUI');
    ud = []; ud.taskState = taskState; ud.beacon_times = beacon_times;
    set(figTC,'UserData',ud);
    % Pause button
    btnPause = uicontrol('Style', 'pushbutton', 'String', getPauseButtonText(), ...
                         'Position', [20, 150, 245, 40], 'Callback', @togglePause, ...
                         'Tag', 'btnPause');

    % Stop Task button
    btnStop = uicontrol('Style', 'pushbutton', 'String', 'Stop Task', ...
                        'Position', [275, 150, 130, 40], 'Callback', @stopTask, ...
                        'Tag', 'btnStop', 'BackgroundColor', [0.85, 0.33, 0.10], ...
                        'ForegroundColor', 'white', 'FontWeight', 'bold');
    
    % Text area for display
    txtDisplay = uicontrol('Style', 'text', 'String', 'Status: Idle', ...
                           'Position', [50, 100, 200, 30], 'BackgroundColor', 'white', ...
                           'FontSize',7,...
                           'Tag', 'txtDisplay');
    
    % Button "Send Event"
    btnSendEvent = uicontrol('Style', 'pushbutton', 'String', 'Send Event', ...
                             'Position', [50, 70, 200, 30], 'Callback', @sendEvent, ...
                             'Tag', 'btnSendEvent', 'BackgroundColor', 'default');
    
    % Button "Send PCPSync"
    btnSendPCPSync = uicontrol('Style', 'pushbutton', 'String', 'Send PCPSync', ...
                               'Position', [50, 40, 200, 30], 'Callback', @sendPCPSync, ...
                               'Tag', 'btnSendPCPSync', 'BackgroundColor', 'default');
    
    % Editable field with text label "Event Code"
    uicontrol('Style', 'text', 'String', 'Event Code:', ...
              'Position', [20, 10, 80, 20], 'BackgroundColor', 'white');
    edtEventCode = uicontrol('Style', 'edit', 'String', '', ...
                             'Position', [110, 10, 140, 20], 'Tag', 'edtEventCode');
    % # [temporary]
    ud.btnPause = btnPause;
    ud.btnStop  = btnStop;
    ud.updateGUIBasedOnTaskState = @ updateGUIBasedOnTaskState;
    set(figTC,'UserData',ud);

    set(btnSendEvent, 'Enable', 'off');
    set(btnStop, 'Enable', 'off');
    set(btnSendPCPSync, 'Enable', 'on');
    set(edtEventCode, 'Enable', 'off');


    function str = getPauseButtonText()
        if ~pause_isActive
            if pause_requested
                str = 'Pause Requsted, press again to cancel ...';
            else
                str = 'Press to Pause';
            end
        else
            if pause_requested
                str = 'Task paused, press again to resume';
            else
                str = 'Resuming...';
            end
        end
    end
    
    function togglePause(~, ~)
        ud = get(figTC,'UserData');
        taskState = ud.taskState;
        pause_requested = taskState.pause_requested;
        %
        pause_requested = ~pause_requested;
        %
        taskState.pause_requested = pause_requested;
        ud.taskState = taskState;
        set(figTC,'UserData',ud)
        updateGUIBasedOnTaskState();
    end
    function updateGUIBasedOnTaskState()
        ud = get(figTC,'UserData'); taskState = ud.taskState;
        pause_requested = taskState.pause_requested;
        pause_isActive  = taskState.pause_isActive;
        stop_requested  = taskState.stop_requested;
        set(btnPause, 'String', getPauseButtonText())
        if pause_requested && ~pause_isActive % requested, waiting for pause
            set(btnPause, 'BackgroundColor', [222,196,159]/255, 'ForegroundColor', 'black');
            set(txtDisplay,'String','<Running> Pause requested, waiting for end of this trial... [Press again to cancel]')
        elseif pause_requested && pause_isActive % task is paused after pause-request
            set(btnPause, 'BackgroundColor', [229,150,150]/255, 'ForegroundColor', 'black');
            set(txtDisplay,'String','<Paused>, press to RESUME \n now you can also press STOP to end the task')
            set(btnStop, 'Enable', 'on');
        elseif ~pause_requested && pause_isActive % task is paused, but resume is requested
            set(btnPause, 'BackgroundColor', [174,235,230]/255, 'ForegroundColor', 'white');
            set(txtDisplay,'String','<Resuming>...')
        else % not paused, no pause requested, default
            set(btnPause, 'BackgroundColor', 'default', 'ForegroundColor', 'black');
            set(txtDisplay,'String','<Running>')
        end
        % --- stop button appearance & confirmation ---
        if stop_requested
            set(btnStop, 'String', 'Stop Requested...', 'BackgroundColor', [0.6, 0.1, 0.0]);
            if pause_isActive
                % Task has paused in response to stop request — prompt user
                answer = questdlg('Stop the task now? This cannot be undone.', ...
                                  'Confirm Stop Task', ...
                                  'Yes, Stop Task', 'Cancel (Stay Paused)', 'Cancel (Stay Paused)');
                if strcmp(answer, 'Yes, Stop Task')
                    stop_task();
                else
                    % Cancelled — clear stop request, leave task paused
                    ud = get(figTC,'UserData');
                    taskState = ud.taskState;
                    taskState.stop_requested = false;
                    taskState.pause_requested = true; % stay paused; user must click resume
                    ud.taskState = taskState;
                    set(figTC,'UserData',ud);
                    set(btnStop, 'String', 'Stop Task', 'BackgroundColor', [0.85, 0.33, 0.10]);
                    set(txtDisplay,'String','<Paused> Stop cancelled — press Resume to continue')
                end
            end
        else
            set(btnStop, 'String', 'Stop Task', 'BackgroundColor', [0.85, 0.33, 0.10]);
        end
    end
    
    function stopTask(~, ~)
        ud = get(figTC,'UserData');
        taskState = ud.taskState;
        taskState.stop_requested = ~taskState.stop_requested;
        if taskState.stop_requested
            % request pause so the task stops at the next safe point
            taskState.pause_requested = true;
        end
        ud.taskState = taskState;
        set(figTC,'UserData',ud);
        updateGUIBasedOnTaskState();
    end

    function sendEvent(~, ~)
        set(btnSendEvent, 'BackgroundColor', 'light yellow');
        drawnow; % Update the UI immediately
        % Simulate the sending process
        send_event('argA', 'argB', 'argC');
        set(btnSendEvent, 'BackgroundColor', 'default');
    end
    
    function sendPCPSync(~, ~)
        set(btnSendPCPSync, 'BackgroundColor', [255, 222, 33]/255);
        drawnow; % Update the UI immediately
        % Simulate the sending process
        beacon_times = [beacon_times, test_Beacon(pulseParam.interval,pulseParam.duration,pulseParam.count)];
        ud.beacon_times = beacon_times;
        set(figTC,'UserData',ud)
        %send_event('argD', 'argE', 'argF'); % Assuming send_event function is reusable
        set(btnSendPCPSync, 'BackgroundColor', 'default');
    end
end

function stop_task()
    % Placeholder — implement task termination logic here
    disp('stop_task() called: terminating task.');
end

function send_event(argA, argB, argC)
    % Placeholder for send_event function
    % Actual implementation would send events to some external system
    disp(['Event sent with arguments: ', argA, ', ', argB, ', ', argC]);
end

function figTC=taskControlGUI(taskState)
    
    pause_requested = taskState.pause_requested;
    pause_isActive = taskState.pause_isActive;
    task_isRunning = taskState.task_isRunning;

    % Initialize global variables if they do not exist
    if isempty(task_isRunning) || ~task_isRunning
        warning('No active FLvoice task running...')
        figTC = [];
        return
    end
    % Create a figure for the GUI
    figTC = figure('Name', 'Task Control', 'NumberTitle', 'off', ...
                 'Position', [100, 100, 400, 200], 'Tag', 'TaskControlGUI');
    ud = []; ud.taskState = taskState;
    set(figTC,'UserData',ud);
    % Pause button
    btnPause = uicontrol('Style', 'pushbutton', 'String', getPauseButtonText(), ...
                         'Position', [50, 150, 200, 40], 'Callback', @togglePause, ...
                         'Tag', 'btnPause');
    
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
    ud.updateGUIBasedOnTaskState = @ updateGUIBasedOnTaskState;
    set(figTC,'UserData',ud);

    set(btnSendEvent, 'Enable', 'off');
    set(btnSendPCPSync, 'Enable', 'off');
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
        pause_isActive = taskState.pause_isActive;
        set(btnPause, 'String', getPauseButtonText())
        if pause_requested && ~pause_isActive % requested, waiting for pause
            set(btnPause, 'BackgroundColor', [222,196,159]/255, 'ForegroundColor', 'black');
            set(txtDisplay,'String','<Running> Pause requested, waiting for end of this trial... [Press again to cancel]')
        elseif pause_requested && pause_isActive % task is paused after pause-request
            set(btnPause, 'BackgroundColor', [229,150,150]/255, 'ForegroundColor', 'black');
            set(txtDisplay,'String','<Paused>, press to RESUME')
        elseif ~pause_requested && pause_isActive % task is paused, but resume is requested
            set(btnPause, 'BackgroundColor', [174,235,230]/255, 'ForegroundColor', 'white');
            set(txtDisplay,'String','<Resuming>...')
        else % not paused, no pause requested, default
            set(btnPause, 'BackgroundColor', 'default', 'ForegroundColor', 'black');
            set(txtDisplay,'String','<Running>')
        end
    end
    
    function sendEvent(~, ~)
        set(btnSendEvent, 'BackgroundColor', 'light yellow');
        drawnow; % Update the UI immediately
        % Simulate the sending process
        send_event('argA', 'argB', 'argC');
        set(btnSendEvent, 'BackgroundColor', 'default');
    end
    
    function sendPCPSync(~, ~)
        set(btnSendPCPSync, 'BackgroundColor', 'light yellow');
        drawnow; % Update the UI immediately
        % Simulate the sending process
        send_event('argD', 'argE', 'argF'); % Assuming send_event function is reusable
        set(btnSendPCPSync, 'BackgroundColor', 'default');
    end
end

function send_event(argA, argB, argC)
    % Placeholder for send_event function
    % Actual implementation would send events to some external system
    disp(['Event sent with arguments: ', argA, ', ', argB, ', ', argC]);
end

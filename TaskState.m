classdef TaskState < handle
    properties (Access = public)
        pause_requested = false;
        pause_isActive = false;
        task_isRunning = false; % "running" is independent of "pause/resume"
        pcpSync_deliveryState = false;
        taskEvent_deliveryState = false;
    end
    events
        StateChanged;  % Event triggered when state changes
        pcpSync_StateChanged; %
        taskEvent_StateCahnged;
    end
    methods
        function set.pause_requested(obj, val)
            obj.pause_requested = val;
            notify(obj, 'StateChanged');
        end
        function set.pause_isActive(obj, val)
            obj.pause_isActive = val;
            notify(obj, 'StateChanged');
        end
        function set.task_isRunning(obj, val)
            obj.task_isRunning = val;
            notify(obj, 'StateChanged');
        end
        function set.pcpSync_deliveryState(obj, val)
            obj.pcpSync_deliveryState = val;
            notify(obj, 'pcpSync_StateChanged');
        end
        function set.taskEvent_deliveryState(obj, val)
            obj.taskEvent_deliveryState = val;
            notify(obj, 'taskEvent_StateCahnged');
        end
    end
end
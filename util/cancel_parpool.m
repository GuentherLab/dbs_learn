function cancel_parpool()

parpool_obj = gcp('nocreate');
if ~isempty(parpool_obj)
    cancel(parpool_obj.FevalQueue.RunningFutures);   % new-style futures queue (R2022b+)
end
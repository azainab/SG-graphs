function [t_events_new] = fEventCount(t_events_old, t_new)

% first check to see if there are "close events" in t_events_old
n = 1; rw_rmv = [];
for qt = 1:length(t_events_old(:,1))-1
    if abs(t_events_old(qt,1)-t_events_old(qt+1,1))<=5
        rw_rmv(n) = qt+1;
        n = n + 1;
    end
end

if ~isempty(rw_rmv)
    t_events_old(rw_rmv,:) = [];
end


n = 1; rw_rmv = [];
for qt = 1:length(t_new)-1
    if abs(t_new(qt)-t_new(qt+1))<=5
        rw_rmv(n) = qt+1;
        n = n + 1;
    end
end

if ~isempty(rw_rmv)
    t_new(rw_rmv) = [];
end

n = 1; new_event =[];
for qt = 1:length(t_new)
    fnd = 0;
    for bb = 1:length(t_events_old(:,1))
        if (t_new(qt) == t_events_old(bb,1))|(t_new(qt) == t_events_old(bb,1)+5)|(t_new(qt) == t_events_old(bb,1)-5)  
            t_events_old(bb,2) = t_events_old(bb,2) + 1;
            fnd = 1;
             break;
        end
    end
    
    if ~fnd
        new_event(n,:) = [t_new(qt), 1];
        n = n + 1;
    end
end

if ~isempty(new_event)
    t_events_new = [t_events_old; new_event];
else
    t_events_new = t_events_old;
end
% using the data from the analysis, find events that are 3 (or other) x standard
% deviation above the average.  t is the time stampe, xAnalysis is the data
% from the analysis, and std_cuttoff is the number of standard deviations
% that is considered an event. the function returns pt_event which
% is a placeholder of where in xAnalysis data the events are  
function [pt_event] = fPlotandAnalysis(timeStamp,signal,t,xAnalysis,std_cuttoff,if_diff,if_plot,sig)


% clear mean_plt std_plt
mean_plt = mean(xAnalysis)*ones(length(xAnalysis),1);
std_plt = std(xAnalysis)*ones(length(xAnalysis),1);

n_pt = 1; pt_event = [];
pt = find(xAnalysis>mean(xAnalysis)+std(xAnalysis)*std_cuttoff);

for qt = 1:length(pt)
    
    t_event = t(pt(qt));
    end_pt = find(t_event==timeStamp);
    if end_pt-10*30==0
        xEvent = signal(end_pt-10*30+1:end_pt);
        tEvent = timeStamp(end_pt-10*30+1:end_pt);
    else
        xEvent = signal(end_pt-10*30:end_pt);
        tEvent = timeStamp(end_pt-10*30:end_pt);
    end
    
    if if_diff
        xEvent_diff = diff(xEvent);
        max_xEvent_diff = max(abs(xEvent_diff));
        pt_max = find(max_xEvent_diff==abs(xEvent_diff));
        t_max = tEvent(pt_max);
    
        if max_xEvent_diff>0.01
            
            if if_plot
                figure;
                plot(tEvent,xEvent), title(['detected event - no. ', num2str(qt), sig])
            end
            
            pt_event(n_pt,1) = pt(qt);
            n_pt = 1 + n_pt;
            
        end
        
    else
        
        if if_plot
            figure;
            plot(tEvent,xEvent), title([' Event no. ',num2str(qt), sig ])
        end
        
        pt_event(n_pt,1) = pt(qt);
        n_pt = 1 + n_pt;
        
    end

end

% plot found events in analyzed data 
if isempty(pt_event)&if_plot
    figure;
    plot(t,xAnalysis,'ro',t,mean_plt,t,std_plt*std_cuttoff+mean_plt,t,mean_plt-std_plt*std_cuttoff), title(['no events', sig])
elseif if_plot
    figure;
    plot(t(pt_event),xAnalysis(pt_event),'kx',t,xAnalysis,'ro',t,mean_plt,t,std_plt*std_cuttoff+mean_plt,t,mean_plt-std_plt*std_cuttoff), title(sig)
end
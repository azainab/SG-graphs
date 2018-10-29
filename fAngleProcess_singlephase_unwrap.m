% Unwraps single phase voltage PMU data
% Input: V1ang - raw PMU phase angle data
% Output: ang - unwrapped PMU phase angle data

function [ang] = fAngleProcess_singlephase_unwrap(V1ang)

angdiff = V1ang(2:end) - V1ang(1:end-1);

posjmp = find(angdiff<-300);
negjmp = find(angdiff>300);

if ~isempty(posjmp)
    for qt = 1:length(posjmp)
        x = [V1ang(1:posjmp(qt)); V1ang(posjmp(qt)+1:end)+360];
%         figure; plot(x)
        V1ang = x;
    end
end

if ~isempty(negjmp)
    for qt = 1:length(negjmp)
        x = [V1ang(1:negjmp(qt)); V1ang(negjmp(qt)+1:end)-360];
%         figure; plot(x)
        V1ang = x;
    end
end

% ang = V1ang - V1ang(1); % start at 0 degrees
ang = V1ang;

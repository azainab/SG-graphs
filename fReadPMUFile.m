% clear all; close all; clc;

%% Read in data and format type
function [timeStamp, stationTitle, stationData, angle_diff_title, angle_difference] = fReadPMUFile(fileName, plts)

[protocol, fileData, fileTitle] = fPMU_protocol(fileName);

%% Number of stations and time stamp (created, not from file)
NoStations = protocol{2,2}(1);
dataLength = length(fileData(:,1));
timeStamp = 0:1/30:1/30*dataLength-1/30;

%% Create Variables for Station Data and Station Data Names
st_rmv = [];
st_error_count = 0;
for qt = 4:6:NoStations*6 %4:6:NoStations*3+4
    st_cnt = 0;
    st = 0;

    for www = 1:3
        qual_test = protocol{2,(qt+(www-1)*2)};
        st(www) = str2double(cell2mat(protocol{2,(qt+(www-1)*2-1)}));
        
        if strcmp(qual_test,'true')
            st_cnt = st_cnt + 1;
        end
        
    end
    
    if st_cnt~=3
%         fileData(:,st(3)+1) = []; fileData(:,st(2)+1) = []; fileData(:,st(1)+1) = [];
%         fileTitle(st(3)+1) = []; fileTitle(st(2)+1) = []; fileTitle(st(1)+1) = [];
%         protocol(:,st(1)*3:st(1)*3+5) = [];
        st_error_count = st_error_count + 1;
        st_rmv = [st_rmv, st];
        st_rmv_protocol(st_error_count) = st(1);
    end
    
end


n = 0;
if ~isempty(st_rmv)
    
    bb = length(st_rmv_protocol);
    aa = length(st_rmv);
    st_rmv_protocol = sort(st_rmv_protocol,'descend');
    st_rmv = sort(st_rmv,'descend');

    for qt = 1:bb
        protocol(:,st_rmv_protocol(qt)*3:st_rmv_protocol(qt)*3+5) = [];
        [anew,bnew] = size(protocol);
        for cnt_proto = st_rmv_protocol(qt)*3:2:bnew-1
            if n<2
                new_loc = str2double(cell2mat(protocol{2,cnt_proto}))-2;
                protocol{2,cnt_proto} = mat2cell(num2str(new_loc));
                n = n + 1;
            else
                new_loc = str2double(cell2mat(protocol{2,cnt_proto}))-3;
                protocol{2,cnt_proto} = mat2cell(num2str(new_loc));
                n = 0;
            end
        end
        
        if st_rmv_protocol(qt)*3~= 3
            for cnt_proto = 7:6:st_rmv_protocol(qt)*3-2
                new_loc = str2double(cell2mat(protocol{2,cnt_proto}))-2;
                protocol{2,cnt_proto} = mat2cell(num2str(new_loc));
            end
        end
    end
    
    for qt = 1:aa
        fileData(:,st_rmv(qt)+1) = []; %fileData(:,st_rmv(qt,2)+1) = []; fileData(:,st_rmv(qt,1)+1) = [];
        fileTitle(st_rmv(qt)+1) = []; %fileTitle(st_rmv(qt,2)+1) = []; fileTitle(st_rmv(qt,1)+1) = [];
    end
end
NoStations = NoStations - st_error_count;

[a,b] = size(protocol);
n = 1;
for qt = 3:2:b
    qual_test = protocol{2,qt+1};

    
    if strcmp(qual_test,'true')
        dataLocation(n,1) = str2double(cell2mat(protocol{2,qt}));
        n = n + 1;
    end
    
end

sv = [];
stationTitle = cell(NoStations*3,1);
for qq = 1:NoStations*3
    temp = fileTitle{dataLocation(qq)+1};
    temp(temp=='_') = ' ';
    stationTitle{qq,1} = temp;
    if strcmp(temp(1:6),'HARRIS')&&isempty(sv)
        sv = qq + 2;
    end
end

stationData = zeros(dataLength,3,NoStations);
for qq = 1:NoStations
    stationData(1:dataLength,1:3,qq) = fileData(:,dataLocation((qq-1)*3+1:qq*3)+1);
end
stationData(1:dataLength,2,sv/3) = stationData(1:dataLength,2,sv/3)+90;

%% Data conditioning
% check if data is "true" first
% [stationData, qual_test] = fdataconditioning(timeStamp, stationData, 0);
% stationTitle(qual_test*3-2:qual_test*3) = [];
% NoStations = NoStations-length(qual_test);


%% Angle Difference Calculation
n = 1;
for bb = NoStations:-1:2
    for qt = 1:bb-1
        angle_diff_title{n,1} = [stationTitle{(qt-1)*3+2,1},' and ',stationTitle{(bb-1)*3+2,1}];
        angle_difference(:,n) = fAngleProcess(stationData(:,2,qt),stationData(:,2,bb));
        
        % TEMP !!!!!!!!!
        rng = 29480:43170;
        N = length(rng);
        
        angle1_unwrap = fAngleProcess_singlephase_unwrap(stationData(rng,2,qt));
        angle2_unwrap = fAngleProcess_singlephase_unwrap(stationData(rng,2,bb));
        
%         figure; plot((1:N)/30,stationData(rng,2,qt),(1:N)/30,stationData(rng,2,bb)), title(angle_diff_title{n,1}),legend(stationTitle{(qt-1)*3+2,1},stationTitle{(bb-1)*3+2,1})
%         figure; plot((1:N)/30,angle1_unwrap,(1:N)/30,angle2_unwrap), title(angle_diff_title{n,1}),legend(stationTitle{(qt-1)*3+2,1},stationTitle{(bb-1)*3+2,1})
%         figure; plot((1:N)/30,angle1_unwrap,(1:N)/30,angle2_unwrap-60), title(angle_diff_title{n,1}),legend(stationTitle{(qt-1)*3+2,1},stationTitle{(bb-1)*3+2,1})
%         figure; plot(angle_difference(rng,n)), title([angle_diff_title{n,1},' difference 1'])
%         figure; plot((1:N)/30,angle2_unwrap-angle1_unwrap), title([angle_diff_title{n,1},' difference 2'])
%         figure; plot((1:N)/30,angle2_unwrap-angle1_unwrap-60), title([angle_diff_title{n,1},' difference 2'])
        
        n = n + 1;
    end
end

%% PLOTS - All Raw Data & Phase Angle Differences
if plts
    n = 1;
    for bb = 1:NoStations
        for qt = 1:3
            figure;
            plot(timeStamp,stationData(:,qt,bb))
            title(stationTitle{n})
            n = n + 1;
        end
    end

    N = sum(1:NoStations-1);
    n = 1;
    for bb = 1:N
            figure;
            plot(timeStamp,angle_difference(:,n)), title(angle_diff_title{n,1})

            n = n + 1;
    end
end
    
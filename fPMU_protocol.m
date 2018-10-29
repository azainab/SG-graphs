function [Protocol, Vmeters, title] = fPMU_protocol(fileName)


% [Vmeters title] = xlsread(fileName);  % Load data

delimiterIn = ',';
headerlinesIn = 1;
A = importdata(fileName,delimiterIn,headerlinesIn);
numdata = A.data;
time_created = (0:1/30:(108000-1)*1/30)';

Vmeters = [time_created, numdata];
title = A.textdata(1,:);
%% TEMP
% if str2num(fileName(8:9))<19
%     Vmeters(:,end-1) = []; title(end-1) = [];
%     Vmeters(:,end-3) = []; title(end-3) = [];
%     Vmeters(:,end-11:end-10) = []; title(end-11:end-10) = [];
% end
%[Vmeters title] = xlsread('Test2.csv');  % Load data


% X = Vmeters(:,1);
% TimeStamp=datevec(X);


%%
% Index of input CSV files. 
% It need to update information if changed or incresed measurement point name
Format_type = {'type1', 'type2','type3'};
index = {'Austin_V1LPM_Magnitude','Austin_V1LPM_Angle','Austin_Frequency','CLOUD_V1LPM_Magnitude','CLOUD_V1LPM_Angle','CLOUD_Frequency','HARRIS_V1LPM_Magnitude','HARRIS_V1LPM_Angle','HARRIS_Frequency','McDonald 1P_V1LPM_Magnitude','McDonald 1P_V1LPM_Angle','McDonald 1P_Frequency','SEL Houston_V1LPM_Magnitude','SEL Houston_V1LPM_Angle','SEL Houston_Frequency','SEL Pullman_V1LPM_Magnitude','SEL Pullman_V1LPM_Angle','SEL Pullman_Frequency','UT Pan Am_V1LPM_Magnitude','UT Pan Am_V1LPM_Angle','UT Pan Am_Frequency'};


[M,N]=size(Vmeters);
Measurement_point = (N-1)/3
Messurement_sample = M
string_cell = cell(2,N-1);

%%
% Find the first frequency value index & define CSV file format
Mean_Vmeters = mean(Vmeters);    % Mean CSV file
Mean_Vmeters = round(Mean_Vmeters);  % floor mean datas
Var_Vmeters = Var(Vmeters);    % Var CSV file
Var_Vmeters = round(Var_Vmeters);  % floor mean datas

for i = 1:1:N
    if Mean_Vmeters(i) == 60
        Position_60Hz(i) = 1;
    else
        Position_60Hz(i) = 0;
    end
 end 
% Mean valuses of frequency are nearly same 60. 

% define CSV file format
% No. 1 -> Mag-Ang-Fre, Mag-Ang-Fre, ...Mag-Ang-Fre
%          0 0 1 | 0 0 1| 0 0 1|......} 0 0 1 
% No. 2 -> Mag-Ang, Mag-Ang, ...Fre-Fre
%          0 0 | 0 0 | 0 0 |......} 1 1 .....  
% No. 3 -> Others

flag_cnt1 = 0;
flag_cnt2 = 0;

% format no. 1 
% accumulate the frequency indesx 
for i = 2:3:N-2
    if Position_60Hz(i:i+2) == [0 0 1]
        flag_cnt1 = flag_cnt1 + 1;
    end    
end  

% format no. 2 
% accumulate the frequency indesx 
for i = 2:2:N-Measurement_point+1
    if Position_60Hz(i:i+1) == [0 0]
        flag_cnt2 = flag_cnt2 + 1;
    end    
end  

% Decision of format type if accumulate number is over 80% about all frequency number.
if flag_cnt1 > Measurement_point*0.8   %% 
    CSV_Format_type = Format_type(1,1);
elseif flag_cnt2 > Measurement_point*0.8   %%         
    CSV_Format_type = Format_type(1,2);
else
    CSV_Format_type = Format_type(1,3);
end 

CSV_Format_type

%%
% check error data, using mean and variance
if strcmp((CSV_Format_type),Format_type(1,1)) == 1     % Case type1 
    error_detector = zeros(1,N-1);
    for i = 2:3:N-2
        temp_Mag = [ Mean_Vmeters(1,i) Var_Vmeters(1,i)];
        temp_Ang = [ Mean_Vmeters(1,i+1) Var_Vmeters(1,i+1)];
        temp_Fre = [ Mean_Vmeters(1,i+2) Var_Vmeters(1,i+2)];
        if temp_Mag == [0 0]
            error_detector(i-1)=1;  
        end
        
        if temp_Ang == [0 0]
            error_detector(i)=1;
        end
        
        if temp_Fre == [0 0]
            error_detector(i+1)=1;
        end
    end
    Number_error_data = sum(error_detector);
    Number_true_data = N-1-Number_error_data;
    
elseif strcmp((CSV_Format_type),Format_type(1,2)) == 1     % Case type2   .
    error_detector = zeros(1,N-1);
    for i = 2:2:N-Measurement_point-1
        temp_Mag = [ Mean_Vmeters(1,i) Var_Vmeters(1,i)];
        temp_Ang = [ Mean_Vmeters(1,i+1) Var_Vmeters(1,i+1)];
        temp_Fre = [ Mean_Vmeters(1,0.5*i+2*Measurement_point+1) Var_Vmeters(1,0.5*i+2*Measurement_point+1)];
  
        if temp_Mag == [0 0]
            error_detector(i-1)=1;
        end
        
        if temp_Ang == [0 0]
            error_detector(i)=1;
        end
        
        if temp_Fre == [0 0]
            error_detector(0.5*i+2*Measurement_point)=1;
        end    
    end    
    Number_error_data = sum(error_detector);
    Number_true_data = N-1-Number_error_data;
end

Number_error_data 
Number_true_data 


%% It just show title and true or falut
for i = 1:1:N-1
    if error_detector(i) == 0
        string_cell(:,i) = [title(i+1); 'true'];
    else
        string_cell(:,i) = [title(i+1); 'fault'];
    end
end


%% Output Protocol
%num2str

if strcmp((CSV_Format_type),Format_type(1,1)) == 1     % Case Type 1
    Protocol = cell(2,2+Measurement_point*6);
    Protocol(:,1) = ['CSV file type'; Format_type(1,1)];
    Protocol(:,2) = ['Number of measurement points'; {Measurement_point} ];
    
    for i = 1:1:Measurement_point
        Protocol(:,6*(i-1)+1+2:6*(i-1)+6+2) = {string_cell(1,3*(i-1)+1), {}, string_cell(1,3*(i-1)+2), {}, string_cell(1,3*(i-1)+3), {} ;...
                                             {num2str(3*(i-1)+1)}, string_cell(2,3*(i-1)+1), {num2str(3*(i-1)+2)}, string_cell(2,3*(i-1)+2), {num2str(3*(i-1)+3)}, string_cell(2,3*(i-1)+3)};   
    end
    
elseif strcmp((CSV_Format_type),Format_type(1,2)) == 1     % Case type2
    Protocol = cell(2,2+Measurement_point*6);
    Protocol(:,1) = ['CSV file type'; Format_type(1,2)];
    Protocol(:,2) = ['Number of measurement points'; {Measurement_point} ];

    for i = 1:1:Measurement_point
        Protocol(:,6*(i-1)+1+2:6*(i-1)+6+2) = {string_cell(1,2*(i-1)+1), {}, string_cell(1,2*(i-1)+2), {}, string_cell(1,i+2*Measurement_point), {} ;...
                                             {num2str(2*(i-1)+1)}, string_cell(2,2*(i-1)+1), {num2str(2*(i-1)+2)}, string_cell(2,3*(i-1)+2), {num2str(i+2*Measurement_point)}, string_cell(2,i+2*Measurement_point)};   
    end 
end


Protocol;



% break 











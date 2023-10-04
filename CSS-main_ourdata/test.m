%% CSS_package (written by Jaehyoung Hong and Hyukpyo Hong)
%% Relieving daytime sleepiness with personalized sleep-wake patterns aligned with the circadian rhythm
clc;
clear;

% User dependant parameter
time_interval = 2;                                                          %Change according to duration of the epoch of users (min), 사용자 기간에 따른 변경(분)
time_start = 24;                                                            %Change according to the starting time of sleep-wake state (h), sleep-wake 상태의 시작 시간에 따른 변화 (h)
tau_c = 24.09;

disp("---CSS_calculation---")
sleep_light = readtable('Input1_sleep_light.csv');                          %Read input1: sleep-wake pattern (1st column) and light exposure (2nd column)
waso_main = readtable('Input2_WASO_main.csv');  

patt = categorical(sleep_light.Var1);                                       %Convert sleep-wake to categorical array
light = table2array(sleep_light(:,2));  

if patt(end) == 'Sleep'
    patt = [patt;'Wake'];
    light = [light;250];
end

waso = table2array(waso_main(:,1));                                         %Get WASO
main_sleep_temp = categorical(waso_main.Var2);                              %Convert main sleep to categorical array
main_sleep = find(main_sleep_temp=='M');

% Parameters
Q_max = 100; theta = 10; sigma = 3; 
mu = 4.2;coef_y = 0.8;const = 1;v_vh = 1.01;coef_x = -0.16;                 %The values of parameters in the mathematical model adopted 
gate = 1;                                                                   %in the computational package

% Baseline simulation before provided sleep-wake pattern to get stable
% sleep-wake pattern entrained to light-dark cycle
start_num = 120;                                                             %Choose 120 days which is usually enough to get stable sleep-wake patterns

if time_start < 12
    time_start = time_start + 24;
end

% Initial light time space 
it_first= linspace(0, 24*start_num, (24*start_num)*30+1);                   %Baseline time_interval is fixed as 2-min
it2 = it_first;
for k = 1 : length(it2)
    it2(k) =  it2(k) - 24 * (floor(it2(k)/24));
end                                                                         %Mod 24 to get exact time

i_first = zeros(1,length(it_first));
for j = 1 : length(i_first)
    if it2(j) >= 6 && it2(j) < 22-0.3
        i_first(j) = 250;
    end
end                                                                         %Light on from 6 a.m. to 21:42 a.m. as 250 lux
                                                                            %Reasonable light on timing according to            
                                                                            %baseline sleep-wake schedule (sleep occurs between 22:00-6:00) 

tspan_first = it_first;                                                     %Time 

V_00 =  [-0.5; -0.25; 5; 0; -11; 14.25];                                    %Initial for simulation (if 6th column is not too big or small, any initial is ok) 
[t0, y0] = ode15s(@(t,V) PCR(t,V,it_first,i_first,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan_first, V_00);
                                                                            %Baseline 120-day simulation
 
                                                                            
y_maximum = max(y0(end-24*30:end,5));                                       
max_index = find(y0(end-24*30:end,5)==y_maximum);                           %Save MA when MA has its highest value for forced wakefulness

y_minimum = min(y0(end-24*30:end,5));
min_index = find(y0(end-24*30:end,5)==y_minimum);                           %Save MA when MA has its lowest value for forced sleep

wake_effort_Vm = y_maximum ;
wake_effort_Vv = min(y0(end-24*30+max_index-1,4)) ;                         %Save VLPO when MA has its highest value for forced wakefulness

sleep_effort_Vm = y_minimum ;
sleep_effort_Vv = max(y0(end-24*30+min_index-1,4)) ;                        %Save VLPO when MA has its lowest value for forced sleep


% Track homeostatic sleep pressure and circadian rhythm according to       
% sleep-wake pattern and light exposure

tspan = zeros(size(patt,1),1);                                              %Base for time of users

tspan(1) = 24*start_num+time_start;                                         %Starting time of users

for j = 2 : length(tspan)
    tspan(j) = tspan(j-1) + time_interval/60;
end                                                                         %Time of users

it = tspan;                                                                 %Light exposure has same timescale with sleep-wake pattern

% Get sleep time / wake time from data 
sleep = find(patt == 'Sleep');                                              %Find sleep_finish time & sleep_num
real_time = find(diff(sleep)~=1);
real_wk_time = sleep(real_time)+1;
real_sl_time = sleep(real_time+1);                                          %Find sleep-starting, finishing time, and #(sleep)


real_wk_time = [real_wk_time;sleep(end)+1];
real_sl_time = [sleep(1);real_sl_time];                                     %Add first sleep-starting, finishing time

real_wk_time_origin = real_wk_time;
real_sl_time_origin = real_sl_time;                                         

% 실제 시간 - waso(수면 후 각성 시작시간)/time_interval으로 변경
for j = 1 : length(real_sl_time)
    if ceil(waso(j,1)/time_interval) ~= 0
    patt(real_wk_time(j,1)-ceil(waso(j,1)/time_interval):real_wk_time(j,1)-1,1)='Wake';
    real_wk_time(j,1) = real_wk_time(j,1) - ceil(waso(j,1)/time_interval);
    end   
end %Decrease wake time as much as waso

for j = 1 : length(real_sl_time)
    real_wk_time(j) = 1 + 8 * (real_wk_time(j) - 1);
    real_sl_time(j) = 1 + 8 * (real_sl_time(j) - 1);
    real_wk_time_origin(j) = 1 + 8 * (real_wk_time_origin(j) - 1);
    real_sl_time_origin(j) = 1 + 8 * (real_sl_time_origin(j) - 1);
end                                                                         %Use more short time interval for model simulation than that of 
                                                                            %actigraphy (Use 1/8)

tspan_temp = zeros(8*(length(tspan)-1) + 1, 1);
tspan_temp(1) = tspan(1);
for j = 1 : length(tspan_temp)
    tspan_temp(j) = tspan(1) + (j-1) * (time_interval/60)/8;
end                                                                         %Use more short time interval for model simulation than that of 
                                                                            %actigraphy (Use 1/8))
tspan = tspan_temp;
st_fi = [1;real_wk_time];

% Between Initial and Day1
ts = 24 * start_num;
it2 = linspace(ts, 24*start_num+time_start, round(8*60/time_interval*(24*start_num+time_start-ts))+1);
it3 = it2;                                                                  %Time between starting time and finshing time of baseline sleep-wake schedule                                                  

 for j = 1 : length(it2)
    if it2(j) >= 24*(start_num)
        it3(j) = it3(j) - 24*start_num;   
    end
end

i2 = zeros(1,length(it3));
for j = 1 : length(i2)
    if it3(j) >= 6
        i2(j) = 250;
    end
end                                                                         %Light on after 6:00

tspan2 = it2;

tspan_total = [tspan_first' ; tspan2(2:end)' ; tspan(2:end)];               %Make total time

i_total = [i_first' ; i2(2:end)' ; light(2:end)];                           %Make total light
it_total = [it_first' ; it2(2:end)' ; it(2:end)];                           %Make time of total light

% Simulation at the gap between the end time of intial 120 days and the
% starting time of actigraphy 
V_0 = y0(end,:)';
[t1_1, y1_1] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan2, V_0);

% When users wear actiwatch very lately we assume no sleep after baseline
% sleep
D_v = -10.2 - (3.37 * 0.5) * ( const + coef_y * y1_1(:,2) + coef_x * y1_1(:,1) ) + v_vh * y1_1(:,6);
wakeup = find(D_v <= 2.46);
if isempty(wakeup) ~= 1
    if isempty(find(D_v(wakeup(1):end) > 2.46, 1)) ~= 1
        sleep_re = find(D_v(wakeup(1):end) > 2.46);
        V_0 =  [y1_1(wakeup(1)+sleep_re(1)-1,1); y1_1(wakeup(1)+sleep_re(1)-1,2); y1_1(wakeup(1)+sleep_re(1)-1,3); wake_effort_Vv; wake_effort_Vm; y1_1(wakeup(1)+sleep_re(1)-1,6)];
        [t1_1(wakeup(1)+sleep_re(1)-1:end,1), y1_1(wakeup(1)+sleep_re(1)-1:end,:)] = ode15s(@(t,V) PCR_shift(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan2(wakeup(1)+sleep_re(1)-1:end), V_0);       
    end
end

flag = 0;                                                                   %Dummy index for checking match between real and predicted sleep-
                                                                            %wake patterns
                                                                            
for day = 1 : length(st_fi)-1
    % Simulation only with light data
    V_0 = y1_1(end,:)';   
    [t1_2,y1_2] = ode15s(@(t,V) PCR(t,V,it_total,i_total,tau_c,mu, v_vh, coef_x, coef_y, const, gate), tspan(st_fi(day):st_fi(day+1),1), V_0);
    a= 0;
while(1)
    a = a + 1;
    matching_temp = Data_matching( day, tau_c, wake_effort_Vm, wake_effort_Vv, sleep_effort_Vm, sleep_effort_Vv,  y1_2, st_fi, patt, t1_2, tspan, real_sl_time, it_total,i_total );
    t1_2 = matching_temp(:,1);
    y1_2 = matching_temp(:,2:7);
    Qm = ( Q_max ./ (1+exp(-(y1_2(:,5)-theta)/sigma)) );
    % Check the difference between data and simulation
    y_temp = zeros(1+(length(y1_2(:,6))-1)/8,1);
    patt_simul = categorical(y_temp(1:end,1),1);
    for j = 1 : length(patt_simul)
        if Qm(1+8*(j-1),1) > 1
            patt_simul(j,1) = 'Wake';                                       %Qm-firing rate of MA population-is closely related with arousal (Qm > 1 : awake)
        elseif Qm(1+8*(j-1),1) <= 1
            patt_simul(j,1) = 'Sleep';                                      %Qm <= 1 : asleep
        end
    end
    diff_patt = find(patt(1+(st_fi(day)-1)/8:1+(st_fi(day+1)-1)/8,1) ~= patt_simul);                  % Compare Wake/Sleep state simulated by model with Wake/Sleep state from actigraphy
    if isempty(diff_patt) == 1 
        break;
    end
    if diff_patt(1) == length(patt_simul)
        break;
    end
    if a == 10
        disp('Too much time, Error occurs');
        flag = 1;
        break
    end
end
    if flag == 1
        break;
    end
    y1_1 = [y1_1(1:end-1,:);y1_2(1:end,:)];
    t1_1 = [t1_1(1:end-1,:);t1_2(1:end,:)];
end

t_total = t1_1(length(tspan2):end,1);
y_total = y1_1(length(tspan2):end,:);
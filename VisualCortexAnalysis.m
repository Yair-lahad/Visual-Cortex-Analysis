clear all
clc
%setting variables
load SpikesX10U12D.mat
units=10;                          % number of neurons
angle_diff=pi/6;                   %difference between angles
last_angle=2*pi-pi/6;              %330
angle_vec=0:angle_diff:last_angle; %vector represent all angles
directions= 12;                    % number of angles tested
repetitions= 200;                  % number of repetitions
bin_len=0.02;                      % bin length
trial_dur=1.28;                    % trial duration
bins=0:bin_len:trial_dur;          % bins vector divided to bins
PSTH=zeros(units,directions,repetitions,length(bins)); %4d array
Spikes_Per_angle=zeros(units,directions,repetitions); % spikes per angle

%% creating Peri-Stimulus Time Histogram (PSTH) data
for neuron=1 : units
   for dir=1 : directions
       for rep=1 : repetitions
           t_vec=SpikesX10U12D(neuron,dir,rep).TimeList; %information of spikes for every rep per angle per neuron
           PSTH(neuron,dir,rep,:)=histc(t_vec,bins);     %adding spikes counts per bin insted of spike times per rep (orgenizing the data difrrently) 
       end
   end   
end

%% Plotting
% Paper parameters
paper_width     = 16.5; %cm
figure_ratio1    = 0.5;
figure_ratio2 =1.1;
% Set font parameters
title_fontsize  = 12;
labels_fontsize = 10;
legend_fontsize = 11;
ticks_fontsize  = 10;
figure ('Color', 'w', 'Units', 'centimeters', 'Position', [1 1 paper_width figure_ratio1*paper_width]); hold on; 
neuron_to_plot=3; % choose which neuron to view
sgtitle('Unit #'+string(neuron_to_plot)+' -PSTH per direction');
for dir=1 : directions
    subplot (2,6,dir); % create subplot on the main figure
    % orgenazing the data - allow us to calculate mean per bin and plot it:
    bins_matrix = reshape(PSTH(neuron_to_plot,dir,:,:), repetitions, length(bins)); 
    rate=mean(bins_matrix); % for each column from bins matrice calculate mean - return vector
    rate=rate/ bin_len;     % dividing by bin duration
    bar(bins,rate,'b');     % creating the bar
    set(gca,'FontSize',ticks_fontsize);
    %setting labels
    if dir>=7
        xlabel('time [sec]','FontSize', labels_fontsize);
    end
    if(dir==1 || dir==7)
        ylabel('rate [HZ]','FontSize', labels_fontsize); 
    end
    ylim([0 25]); xlim([0 1.28]);            %setting limits for subplot
    title(strcat('\theta','='," ",num2str(rad2deg((dir-1)*angle_diff)),'\circ'),'FontSize', title_fontsize) %representing each angle
end
hold off;


%% 2  Orientation and direction tuning
Spikes_Per_angle=sum(PSTH,4);                  %sum of all spikes of all repetitions of each unit and angle
Rate_Neuron_angle= Spikes_Per_angle/trial_dur; %calculating firing rate for each neuron and angle
Mean_rate= mean(Rate_Neuron_angle,3);          %calculating the mean of the firing rate
Std_rate= std(Rate_Neuron_angle,0,3);          %calculating the std of the firing rate
[Po,idx]=max(Mean_rate');                      %maximum firing rate per preffered oreintation and thier indexes
fitResult=cell(10,1);                          %saving the best option of fit per neuron
X_vec = 0:0.01:2*pi;                           %values used for fit plotting
funType=zeros(1,10);                           %1 represents VM_Direction, 0 repesents VM_Orentation
% direction selective:
% po- preffered orientation
 VM_drct = 'A * exp (k * cos (x - PO))'; % A, k & PO are fitted for independent variable x
% orientation selective:
 VM_ornt = 'A * exp (k * cos (2*(x- PO)))';
% Define the independent variable and the coefficients to be fitted
FitDeff_D = fittype(VM_drct, ...
                  'coefficients', {'A', 'PO','k'}, ...
                  'independent', 'x');
FitDeff_O = fittype(VM_ornt, ...
                  'coefficients', {'A', 'PO','k'}, ...
                  'independent', 'x');
% Define the coefficients' exploration space
A_upper=max(Po); %setting upper range for coefficient A
fitOpt_D = fitoptions (FitDeff_D);
fitOpt_D.Lower       = [0     , 0 , 0];
fitOpt_D.Upper       = [A_upper , pi*2,inf];

fitOpt_O = fitoptions (FitDeff_O);
fitOpt_O.Lower       =  [0     , 0 ,0];
fitOpt_O.Upper       = [A_upper , pi*2,inf];

%% calculating Fit per neuron
for neuron=1:units
    %setting startpoint for each function: 
    %Po-maximum firing rate at Preffered orientation, idx- the angle Po represents
    fitOpt_D.Startpoint  = [Po(neuron),(idx(neuron)-1)*angle_diff,0.5]; 
    fitOpt_O.Startpoint  = [Po(neuron),(idx(neuron)-1)*angle_diff,0.5];
    [fitResultD,GoF_D]= fit(angle_vec',Mean_rate(neuron,:)',FitDeff_D,fitOpt_D);     % GOF = goodness of fit for direction
    [fitResultO, GoF_O] = fit(angle_vec' ,Mean_rate(neuron,:)', FitDeff_O,fitOpt_O); % GOf for orientaion
    if(GoF_D.rmse<GoF_O.rmse) % if rmse is smaller, fit works better
        fitResult{neuron}=fitResultD;
        funType(neuron)=1;
    else
        fitResult{neuron}=fitResultO;
        
    end      
end

%% plotting
figure ('Color', 'w', 'Units', 'centimeters', 'Position', [0 0 paper_width figure_ratio2*paper_width]); hold on;
sgtitle('Direction/orientation selectivity - von Mises fit per unit');
x_ticks=0:90:360; % representing x
for unit_idx = 1:10
    subplot (4,3,unit_idx); hold on; % create subplot on the main figure
    % plot the experimental data:
    errorbar (rad2deg(angle_vec), Mean_rate(unit_idx,:), Std_rate(unit_idx,:), 'o');
    % plot the fitted curve :
    plot (rad2deg(X_vec),fitResult{unit_idx}(X_vec), 'r');
    ax = gca;              % Get currect axes
    ax.FontSize = ticks_fontsize;
    %setting labels
    if unit_idx>7
        xlabel('direction [deg]','FontSize', labels_fontsize);
    end
    if(mod(unit_idx,3)==1)
        ylabel('rate [HZ]','FontSize', labels_fontsize);
    end
    xticks(x_ticks);
    xlim([0 360]);
    title(strcat('Unit #',num2str(unit_idx))) %representing each unit
end
legend('rate','VM fit','FontSize',legend_fontsize,'Position',[0.7 0.125 0.2 0.1]);
    

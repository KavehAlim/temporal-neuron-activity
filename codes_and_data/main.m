%% load data
load('spikes.mat');

%% construct large spikes vector and large spikes matrix
max_sample_num = 0;
for i = 1:size(block.segments{1,1}.spiketrains,2)
   if(max( block.segments{1,1}.spiketrains{1,i}.times ) >= max_sample_num )
       max_sample_num = max( block.segments{1,1}.spiketrains{1,i}.times );
   end
end

%construct a large spikes vector that contains all spike events
large_spikes_vector = zeros(1, max_sample_num);

for i = 1:size(block.segments{1,1}.spiketrains,2)
   for j = 1:size(block.segments{1,1}.spiketrains{1,i}.times, 2)
       large_spikes_vector(block.segments{1,1}.spiketrains{1,i}.times(j)) = 1;
   end
end

save('large_spikes_vector', 'large_spikes_vector');

%% extract trial times
clc;
event_labels = block.segments{1,1}.events{1,1}.labels;
event_times = block.segments{1,1}.events{1,1}.times;

trial_start_finish = [];
for i = 1:length(event_times)-9
   if( str2double(event_labels(i, :)) == 65296 )
       trial_start_finish = [trial_start_finish, event_times(i), event_times(i+9)];
   end
end

%% construct raster matrix
clc;
load('large_spikes_vector.mat');

trial_max_duration = 0;
for i = 1:2:length(trial_start_finish)-1
   if( trial_start_finish(i+1) - trial_start_finish(i) >= trial_max_duration )
      trial_max_duration =  trial_start_finish(i+1) - trial_start_finish(i);
   end
end

trial_spikes = zeros(length(trial_start_finish)/2, trial_max_duration);
for i = 1:2:length(trial_start_finish)
   trial_spikes((i-1)/2+1, 1:trial_start_finish(i+1)-trial_start_finish(i)) = large_spikes_vector(trial_start_finish(i):trial_start_finish(i+1)-1); 
end

%% calculate average SR-onset times
SRONSET_avg = 0;
for i = 1:length(event_times)
    if(str2double(event_labels(i, :)) == 65382 || str2double(event_labels(i, :)) == 65385)
       SRONSET_avg = SRONSET_avg + (event_times(i)-event_times(i-5))/161; 
    end
end

%% construct time-locked spikes data
fs = 30000;
binsize = 1000;
time_locked_vector_length = (1.5 + 0.8) * fs;
time_locked_matrix = [];

for i = 1:length(event_times)
    if(str2double(event_labels(i, :)) == 65366 || str2double(event_labels(i, :)) == 65369)
       time_locked_matrix = [time_locked_matrix; large_spikes_vector(event_times(i) - 0.8*fs : event_times(i) + 1.5*fs-1)];
    end
end

%% calculate and plot time-locked PSTH
clc;
time_locked_PSTH_matrix = moving_average_m(time_locked_matrix, binsize);

time_locked_PSTH = zeros(1, size(time_locked_PSTH_matrix,2));

for i = 1:size(time_locked_PSTH_matrix,1)
    time_locked_PSTH = time_locked_PSTH + (1/size(time_locked_PSTH_matrix,1))*time_locked_PSTH_matrix(i,:);
end

figure;
    plot(linspace(-0.8, 1.5, length(time_locked_PSTH)), time_locked_PSTH * fs, 'k');
    title('PSTH of time-locked spikes'); xlabel('time relative to the onset(s)'); xlim([-0.8 1.5]);
    xline(0, 'r-', 'LineWidth', 2); hold on;
    legend('time-locked PSTH curve', 'GO-ON'); ylabel('cumulative firing rate');
    
%% plot general raster plots
clc;
fs = 30000;
downsample_rate = 100;
% PGLF
trial_spikes_donsampled = (downsample((trial_spikes.'), downsample_rate)).' ;
figure;
for i = 1:size(trial_spikes,1)
    scatter(linspace(0, size(trial_spikes,2)/fs, size(trial_spikes_donsampled,2)), (size(trial_spikes,1)-i+1)*trial_spikes_donsampled(i, :), 'k.'); hold on;
    title('all trials spikes raster plot'); xlabel('task time(s)');
end
    xline(0.8, 'r-', 'LineWidth', 2); hold on; % CUE-ON
    xline(2.1, 'b-', 'LineWidth', 2); hold on; % GO-ON
    xline(0.4, 'c-', 'LineWidth', 2); hold on; % WS-ON
    xline(SRONSET_avg/fs, 'g-', 'Linewidth', 2); % SR-ON
    ylabel('trial number');

%% general PSTH

trial_min_duration = trial_start_finish(2) - trial_start_finish(1);
for i = 1:2:length(trial_start_finish)-1
   if( trial_start_finish(i+1) - trial_start_finish(i) <= trial_min_duration )
      trial_min_duration =  trial_start_finish(i+1) - trial_start_finish(i);
   end
end

binsize = 1000;
trial_PSTH_matrix = moving_average_m(trial_spikes, binsize);

trial_PSTH = zeros(1, size(trial_PSTH_matrix,2));

for i = 1:size(trial_PSTH_matrix,1)
    trial_PSTH = trial_PSTH + (1/size(trial_PSTH_matrix,1))*trial_PSTH_matrix(i,:);
end

figure;
    plot(linspace(0, length(trial_PSTH)/fs, length(trial_PSTH)), trial_PSTH*fs, 'k');
    title('PSTH of all trial spikes'); xlabel('task time(s)'); xlim([0,length(trial_PSTH)/fs]);
    xline(0.8, 'r-', 'LineWidth', 1.5); hold on; 
    xline(2.1, 'b-', 'LineWidth', 1.5); hold on;
    xline(0.4, 'c-', 'LineWidth', 1.5); hold on;
    xline(SRONSET_avg/fs, 'g-', 'LineWidth', 1.5); xlim([0, trial_min_duration/fs]); legend({'PSTH curve', 'CUE-ON', 'GO-ON', 'WS-ON', 'SR-ON'}, 'Location', 'southeast'); ylabel('cumulative firing rate (Hz)');
    
%% most eventful events!
clc;
rates_before_GO = [];
rates_after_GO = [];
rates_before_SR = [];
rates_after_SR = [];

averaging_binsize = 1 * fs;
for i = 1:length(event_times)
    
    if( str2double(event_labels(i, :)) == 65366 || str2double(event_labels(i, :)) == 65369 ) % rates before and after GO-ON
       before_avg = mean( large_spikes_vector(event_times(i)-averaging_binsize:event_times(i)) );
       after_avg = mean( large_spikes_vector(event_times(i):event_times(i)+averaging_binsize) );
       rates_before_GO = [rates_before_GO, before_avg];
       rates_after_GO = [rates_after_GO, after_avg];
    end
    
    if( str2double(event_labels(i, :)) == 65382 || str2double(event_labels(i, :)) == 65385 ) % rates before and after SR-ON
       before_avg = mean( large_spikes_vector(event_times(i)-averaging_binsize:event_times(i)) );
       after_avg = mean( large_spikes_vector(event_times(i):event_times(i)+averaging_binsize) );
       rates_before_SR = [rates_before_SR, before_avg];
       rates_after_SR = [rates_after_SR, after_avg];
    end  
end

SR_rate_difference = rates_after_SR - rates_before_SR;
GO_rate_difference = rates_after_GO - rates_before_GO;

figure;
    subplot(211);
        histogram(SR_rate_difference*fs, floor(length(SR_rate_difference)/2), 'FaceColor', [0.5 0 0.1]); title('histogram of rate difference after SR-ON'); xlabel('cumulative firing rate difference (Hz)'); ylabel('population');
        xlim([0, 950]);
    subplot(212);
        histogram(GO_rate_difference*fs, floor(length(GO_rate_difference)/2), 'FaceColor', [0.1 0.5 0]); title('histogram of rate difference after GO-ON'); xlabel('cumulative firing rate difference (Hz)'); ylabel('population');
        xlim([0, 950]);
        
%% apply t-test
% for GO-ON set
    mu = 0; % expected change in firing rate for null hypothesis
    [h_GO, p_GO] = ttest(SR_rate_difference-mu);
    GO_t_val = (mean(SR_rate_difference)-mu) / (std(SR_rate_difference, 0, 'all')/sqrt(length(SR_rate_difference)));

    % plot the t distribution
        v = length(SR_rate_difference) - 1; % degrees of freedom
        t = linspace(-70, 70, 7000);
        df = ( gamma((v+1)/2)*(1+t.^2/v).^(-(v+1)/2) ) / (sqrt(v*pi)*gamma(v/2));
        figure;
        subplot(211);
            plot(t, df, 'LineWidth', 1);
            xline(GO_t_val, 'r', 'LineWidth', 1);
            title('t-distribution PDF for GO-ON set');
            xlabel('t');
            legend('PDF', 't-value');
        
% for SR-ON set
    mu = 0; % expected change in firing rate for null hypothesis
    [h_SR, p_SR] = ttest(GO_rate_difference-mu);
    SR_t_val = (mean(GO_rate_difference)-mu) / (std(GO_rate_difference, 0, 'all')/sqrt(length(GO_rate_difference)));
    
    % plot the t distribution
        v = length(GO_rate_difference) - 1; % degrees of freedom
        t = linspace(-70, 70, 7000);
        df = ( gamma((v+1)/2)*(1+t.^2/v).^(-(v+1)/2) ) / (sqrt(v*pi)*gamma(v/2));
        subplot(212);
            plot(t, df, 'LineWidth', 1);
            xline(SR_t_val, 'r', 'LineWidth', 1);
            title('t-distribution PDF for SR-ON set');
            xlabel('t');
            legend('PDF', 't-value');
            
%% fano-factor and change in distribution
clc;
fano_before_event = [];
fano_after_event = [];
GO_times = [];

for i = 1:length(event_times)
   if( str2double(event_labels(i, :)) == 65366 || str2double(event_labels(i, :)) == 65369 )
      GO_times = [GO_times, event_times(i)];
   end
end

max_trial_duration = 0;
for i = 1:length(block.segments{1,1}.spiketrains)
    if(block.segments{1,1}.spiketrains{1,i}.times(length(block.segments{1,1}.spiketrains{1,i}.times)) >= max_trial_duration)
       max_trial_duration =  block.segments{1,1}.spiketrains{1,i}.times(length(block.segments{1,1}.spiketrains{1,i}.times));
    end
end

for i = 1:length(GO_times) % per trial
   for j = 1: length(block.segments{1,1}.spiketrains) % per neuron
      
      neuron_spikes = zeros(1, max_trial_duration);
      
      neuron_spikes(block.segments{1,1}.spiketrains{1,j}.times) = 1;
      
      firing_rate_before = mean(neuron_spikes(GO_times(i)-0.1*fs : GO_times(i)));
      firing_rate_after = mean(neuron_spikes(GO_times(i) : GO_times(i)+0.1*fs));
      
      psth_before = moving_average_m(neuron_spikes(GO_times(i)-0.1*fs : GO_times(i)), binsize);
      psth_after = moving_average_m(neuron_spikes(GO_times(i) : GO_times(i)+0.1*fs), binsize);
      
      if((firing_rate_after) - (firing_rate_before) > 0 && mean(psth_after)~=0 && mean(psth_before) ~= 0)
          fano_after_event =[fano_after_event, var(psth_after)/mean(psth_after)];
          fano_before_event = [fano_before_event, var(psth_before)/mean(psth_before)];
      end
   end
end
save('rate_fano_after', 'fano_after_event');
save('rate_fano_before', 'fano_before_event');
%% plot fano factor recordings
load('rate_fano_after.mat');
load('rate_fano_before.mat');

fano_after_std = sqrt(var(fano_after_event(~isnan(fano_after_event))));
fano_after_mean = mean(fano_after_event(~isnan(fano_after_event)));
fano_before_std = sqrt(var(fano_before_event(~isnan(fano_before_event))));
fano_before_mean = mean(fano_before_event(~isnan(fano_before_event)));

figure;
    histogram(fano_after_event(~isnan(fano_after_event)), floor(length(fano_after_event)/100)); hold on;
    histogram(fano_before_event(~isnan(fano_before_event)), floor(length(fano_before_event)/100)); 
    legend('after GO-ON event', 'before GO-ON event');
    title('comparisson of fano factor distribution before and after GO-ON event'); xlabel('fano factor (f)'); ylabel('population');
    
%% apply t-test on fano factor data
clc;
mu = 0; % expected change in firing rate for null hypothesis
fano_difference = fano_after_event - fano_before_event;
fano_difference = fano_difference(~isnan(fano_difference));
% fano_difference_std = sqrt(var(fano_difference));
% fano_difference_mean = mean(fano_difference);
% fano_difference = fano_difference(abs(fano_difference-fano_difference_mean) <= 2*fano_difference_std);
[h_fano, p_fano] = ttest(fano_difference-mu);
fano_t_val = (mean(fano_difference)-mu) / (std(fano_difference, 0, 'all')/sqrt(length(fano_difference)));

% plot the t distribution
    v = length(fano_difference) - 1; % degrees of freedom
    t = linspace(-30, 30, 3000);
    df = ( gamma((v+1)/2)*(1+t.^2/v).^(-(v+1)/2) ) / (sqrt(v*pi)*gamma(v/2));
    figure;
        plot(t, df, 'LineWidth', 1);
%         xline(fano_t_val, 'r', 'LineWidth', 1);
        title('t-distribution PDF for fano factor recordings');
        xlabel('t');
        legend('PDF', 't-value');

%% plot fano difference distribution
fano_difference_std = sqrt(var(fano_difference));
fano_difference_mean = mean(fano_difference);
histogram(fano_difference(abs(fano_difference-fano_difference_mean) <= 3*fano_difference_std), floor(length(fano_difference)/100));
title('histogram of fano factor changes after GO onset'); xlabel('change in fano factor'); ylabel('population');

%% extract trial types and their respective times

PGLF_times = [];
PGHF_times = [];
SGLF_times = [];
SGHF_times = [];

for i = 1:length(event_times)
    % PG HF times [TS-ON SR STOP]
    if(str2double(event_labels(i, :)) == 65365 && str2double(event_labels(i+2, :)) == 65366)
       PGHF_times = [PGHF_times, event_times(i-2), event_times(i+3), event_times(i+6)]; 
    end
    
    % PG LF times [TS-ON SR STOP]
    if(str2double(event_labels(i, :)) == 65365 && str2double(event_labels(i+2, :)) == 65369)
       PGLF_times = [PGLF_times, event_times(i-2), event_times(i+3), event_times(i+6)]; 
    end
    
    % SG HF times [TS-ON SR STOP]
    if(str2double(event_labels(i, :)) == 65370 && str2double(event_labels(i+2, :)) == 65366)
       SGHF_times = [SGHF_times, event_times(i-2), event_times(i+3), event_times(i+6)]; 
    end
    
    % SG LF times [TS-ON SR STOP]
    if(str2double(event_labels(i, :)) == 65370 && str2double(event_labels(i+2, :)) == 65369)
       SGLF_times = [SGLF_times, event_times(i-2), event_times(i+3), event_times(i+6)]; 
    end
    
end

save('PGLF_times', 'PGLF_times');
save('PGHF_times', 'PGHF_times');
save('SGLF_times', 'SGLF_times');
save('SGHF_times', 'SGHF_times');

%% categorized raster plot matrices
clc;

% calculate matrix size needed for event spikes recording
PGLF_max_duration = 0;
PGHF_max_duration = 0;
SGLF_max_duration = 0;
SGHF_max_duration = 0;
for i = 1:3:length(PGLF_times)
   if(PGLF_times(i+2)-PGLF_times(i) >= PGLF_max_duration)
       PGLF_max_duration = PGLF_times(i+2)-PGLF_times(i);
   end
end
for i = 1:3:length(PGHF_times)
   if(PGHF_times(i+2)-PGHF_times(i) >= PGHF_max_duration)
       PGHF_max_duration = PGHF_times(i+2)-PGHF_times(i);
   end
end
for i = 1:3:length(SGLF_times)
   if(SGLF_times(i+2)-SGLF_times(i) >= SGLF_max_duration)
       SGLF_max_duration = SGLF_times(i+2)-SGLF_times(i);
   end
end
for i = 1:3:length(SGHF_times)
   if(SGHF_times(i+2)-SGHF_times(i) >= SGHF_max_duration)
       SGHF_max_duration = SGHF_times(i+2)-SGHF_times(i);
   end
end

% record spikes related to each event type
PGLF_spikes = zeros(length(PGLF_times)/3, PGLF_max_duration);
PGHF_spikes = zeros(length(PGHF_times)/3, PGHF_max_duration);
SGLF_spikes = zeros(length(SGLF_times)/3, SGLF_max_duration);
SGHF_spikes = zeros(length(SGHF_times)/3, SGHF_max_duration);
for i = 1:3:length(PGLF_times)
   PGLF_spikes((i-1)/3+1, 1:PGLF_times(i+2)-PGLF_times(i)) = large_spikes_vector(PGLF_times(i):PGLF_times(i+2)-1); 
end
for i = 1:3:length(PGHF_times)
   PGHF_spikes((i-1)/3+1, 1:PGHF_times(i+2)-PGHF_times(i)) = large_spikes_vector(PGHF_times(i):PGHF_times(i+2)-1); 
end
for i = 1:3:length(SGLF_times)
   SGLF_spikes((i-1)/3+1, 1:SGLF_times(i+2)-SGLF_times(i)) = large_spikes_vector(SGLF_times(i):SGLF_times(i+2)-1); 
end
for i = 1:3:length(SGHF_times)
   SGHF_spikes((i-1)/3+1, 1:SGHF_times(i+2)-SGHF_times(i)) = large_spikes_vector(SGHF_times(i):SGHF_times(i+2)-1); 
end


%% plot categorized raster plots
clc;
fs = 30000;
downsample_rate = 100;
% PGLF
PGLF_spikes_donsampled = (downsample((PGLF_spikes.'), downsample_rate)).' ;
figure;
for i = 1:size(PGLF_spikes,1)
    scatter(linspace(0, size(PGLF_spikes,2)/fs, size(PGLF_spikes_donsampled,2)), (size(PGLF_spikes,1)-i+1)*PGLF_spikes_donsampled(i, :), 'k.'); hold on;
    title('PGLF spikes raster plot'); xlabel('t(s)');
end
    xline(0.8, 'r-', 'LineWidth', 2); hold on; 
    xline(2.1, 'b-', 'LineWidth', 2); hold on;
    xline(0.4, 'c-', 'LineWidth', 2); hold on;
    xline(SRONSET_avg/fs, 'g-', 'Linewidth', 2); ylabel('trial number'); xlabel('task time(s)');

% PGHF
PGHF_spikes_donsampled = (downsample((PGHF_spikes.'), downsample_rate)).' ;
figure;
for i = 1:size(PGHF_spikes,1)
    scatter(linspace(0, size(PGHF_spikes,2)/fs, size(PGHF_spikes_donsampled,2)), (size(PGHF_spikes,1)-i+1)*PGHF_spikes_donsampled(i, :), 'k.'); hold on;
    title('PGHF spikes raster plot'); xlabel('t(s)');
end
    xline(0.8, 'r-', 'LineWidth', 2); hold on; 
    xline(2.1, 'b-', 'LineWidth', 2); hold on;
    xline(0.4, 'c-', 'LineWidth', 2); hold on;
    xline(SRONSET_avg/fs, 'g-', 'Linewidth', 2); ylabel('trial number'); xlabel('task time(s)');
    
% SGLF
SGLF_spikes_donsampled = (downsample((SGLF_spikes.'), downsample_rate)).' ;
figure;
for i = 1:size(SGLF_spikes,1)
    scatter(linspace(0, size(SGLF_spikes,2)/fs, size(SGLF_spikes_donsampled,2)), (size(SGLF_spikes,1)-i+1)*SGLF_spikes_donsampled(i, :), 'k.'); hold on;
    title('SGLF spikes raster plot'); xlabel('t(s)');
end
    xline(0.8, 'r-', 'LineWidth', 2); hold on; 
    xline(2.1, 'b-', 'LineWidth', 2); hold on;
    xline(0.4, 'c-', 'LineWidth', 2); hold on;
    xline(SRONSET_avg/fs, 'g-', 'Linewidth', 2); ylabel('trial number'); xlabel('task time(s)');

% SGHF
SGHF_spikes_donsampled = (downsample((SGHF_spikes.'), downsample_rate)).' ;
figure;
for i = 1:size(SGHF_spikes,1)
    scatter(linspace(0, size(SGHF_spikes,2)/fs, size(SGHF_spikes_donsampled,2)), (size(SGHF_spikes,1)-i+1)*SGHF_spikes_donsampled(i, :), 'k.'); hold on;
    title('SGHF spikes raster plot'); xlabel('t(s)');
end
    xline(0.8, 'r-', 'LineWidth', 2); hold on; 
    xline(2.1, 'b-', 'LineWidth', 2); hold on;
    xline(0.4, 'c-', 'LineWidth', 2); hold on;
    xline(SRONSET_avg/fs, 'g-', 'Linewidth', 2); ylabel('trial number'); xlabel('task time(s)');

%% categorized PSTH plots
binsize = 1000;
PGLF_PSTH_matrix = moving_average_m(PGLF_spikes, binsize);
PGHF_PSTH_matrix = moving_average_m(PGHF_spikes, binsize);
SGLF_PSTH_matrix = moving_average_m(SGLF_spikes, binsize);
SGHF_PSTH_matrix = moving_average_m(SGHF_spikes, binsize);

PGLF_PSTH = zeros(1, size(PGLF_PSTH_matrix,2));
PGHF_PSTH = zeros(1, size(PGHF_PSTH_matrix,2));
SGLF_PSTH = zeros(1, size(SGLF_PSTH_matrix,2));
SGHF_PSTH = zeros(1, size(SGHF_PSTH_matrix,2));

for i = 1:size(PGLF_PSTH_matrix,1)
    PGLF_PSTH = PGLF_PSTH + (1/size(PGLF_PSTH_matrix,1))*PGLF_PSTH_matrix(i,:);
end
for i = 1:size(PGHF_PSTH_matrix,1)
    PGHF_PSTH = PGHF_PSTH + (1/size(PGHF_PSTH_matrix,1))*PGHF_PSTH_matrix(i,:);
end
for i = 1:size(SGLF_PSTH_matrix,1)
    SGLF_PSTH = SGLF_PSTH + (1/size(SGLF_PSTH_matrix,1))*SGLF_PSTH_matrix(i,:);
end
for i = 1:size(SGHF_PSTH_matrix,1)
    SGHF_PSTH = SGHF_PSTH + (1/size(SGHF_PSTH_matrix,1))*SGHF_PSTH_matrix(i,:);
end

%% plot PSTH of different events

figure;
    subplot(221);
        plot(linspace(0, length(PGLF_PSTH)/fs, length(PGLF_PSTH)), PGLF_PSTH*fs, 'k');
        title('PSTH of PGLF spikes'); xlabel('task time(s)'); ylabel('cumulative firing rate (Hz)');
        xline(0.8, 'r-', 'LineWidth', 1.5); hold on; 
        xline(2.1, 'b-', 'LineWidth', 1.5); hold on;
        xline(0.4, 'c-', 'LineWidth', 1.5); hold on;
        xline(SRONSET_avg/fs, 'g-', 'LineWidth', 1.5);xlim([0, trial_min_duration/fs]); legend({'PSTH curve', 'CUE-ON', 'GO-ON', 'WS-ON', 'SR-ON'}, 'Location', 'southeast'); 
        ylim([0, 0.12]*fs);
    subplot(222);
        plot(linspace(0, length(PGHF_PSTH)/fs, length(PGHF_PSTH)), PGHF_PSTH*fs, 'k');
        title('PSTH of PGHF spikes'); xlabel('task time(s)'); ylabel('cumulative firing rate (Hz)');
        xline(0.8, 'r-', 'LineWidth', 1.5); hold on; 
        xline(2.1, 'b-', 'LineWidth', 1.5); hold on;
        xline(0.4, 'c-', 'LineWidth', 1.5); hold on;
        xline(SRONSET_avg/fs, 'g-', 'LineWidth', 1.5);xlim([0, trial_min_duration/fs]); legend({'PSTH curve', 'CUE-ON', 'GO-ON', 'WS-ON', 'SR-ON'}, 'Location', 'southeast');
        ylim([0, 0.12]*fs);
    subplot(223);
        plot(linspace(0, length(SGLF_PSTH)/fs, length(SGLF_PSTH)), SGLF_PSTH*fs, 'k');
        title('PSTH of SGLF spikes'); xlabel('task time(s)'); ylabel('cumulative firing rate (Hz)');
        xline(0.8, 'r-', 'LineWidth', 1.5); hold on; 
        xline(2.1, 'b-', 'LineWidth', 1.5); hold on;
        xline(0.4, 'c-', 'LineWidth', 1.5); hold on;
        xline(SRONSET_avg/fs, 'g-', 'LineWidth', 1.5);xlim([0, trial_min_duration/fs]); legend({'PSTH curve', 'CUE-ON', 'GO-ON', 'WS-ON', 'SR-ON'}, 'Location', 'southeast');
        ylim([0, 0.12]*fs);
    subplot(224);
        plot(linspace(0, length(SGHF_PSTH)/fs, length(SGHF_PSTH)), SGHF_PSTH*fs, 'k');
        title('PSTH of SGHF spikes'); xlabel('task time(s)'); ylabel('cumulative firing rate (Hz)');
        xline(0.8, 'r-', 'LineWidth', 1.5); hold on; 
        xline(2.1, 'b-', 'LineWidth', 1.5); hold on;
        xline(0.4, 'c-', 'LineWidth', 1.5); hold on;
        xline(SRONSET_avg/fs, 'g-', 'LineWidth', 1.5);xlim([0, trial_min_duration/fs]); legend({'PSTH curve', 'CUE-ON', 'GO-ON', 'WS-ON', 'SR-ON'}, 'Location', 'southeast');
        ylim([0, 0.12]*fs);
        
%% in one plot for comparisson

plot(linspace(0, length(PGLF_PSTH)/fs, length(PGLF_PSTH)), PGLF_PSTH*fs); hold on;
plot(linspace(0, length(PGHF_PSTH)/fs, length(PGHF_PSTH)), PGHF_PSTH*fs); hold on;
plot(linspace(0, length(SGLF_PSTH)/fs, length(SGLF_PSTH)), SGLF_PSTH*fs); hold on;
plot(linspace(0, length(SGHF_PSTH)/fs, length(SGHF_PSTH)), SGHF_PSTH*fs); hold on;

        title('PSTH of event-based spikes'); xlabel('task time(s)');
        xline(0.8, 'r-', 'LineWidth', 2); hold on; 
        xline(2.1, 'b-', 'LineWidth', 2); hold on;
        xline(0.4, 'c-', 'LineWidth', 2); hold on;
        xline(SRONSET_avg/fs, 'g-', 'LineWidth', 2);xlim([0, trial_min_duration/fs]); legend({'PGLF PSTH', 'PGHF PSTH', 'SGLF PSTH', 'SGHF PSTH', 'CUE-ON', 'GO-ON', 'WS-ON', 'SR-ON'}, 'Location', 'southeast'); ylabel('cumulative firing rate (Hz)');
        
%% PSTH on single neurons
clc;
neuron_id = randi([1, length(block.segments{1,1}.spiketrains)]);
% neuron_id = 80;
spikes = zeros(size(large_spikes_vector));
for i = 1:length(block.segments{1,1}.spiketrains{1,neuron_id}.times)
   spikes( block.segments{1,1}.spiketrains{1,neuron_id}.times(i) ) = 1;
end

% record spikes related to each event type
single_PGLF_spikes = zeros(length(PGLF_times)/3, PGLF_max_duration);
single_PGHF_spikes = zeros(length(PGHF_times)/3, PGHF_max_duration);
single_SGLF_spikes = zeros(length(SGLF_times)/3, SGLF_max_duration);
single_SGHF_spikes = zeros(length(SGHF_times)/3, SGHF_max_duration);
for i = 1:3:length(PGLF_times)
   single_PGLF_spikes((i-1)/3+1, 1:PGLF_times(i+2)-PGLF_times(i)) = spikes(PGLF_times(i):PGLF_times(i+2)-1); 
end
for i = 1:3:length(PGHF_times)
   single_PGHF_spikes((i-1)/3+1, 1:PGHF_times(i+2)-PGHF_times(i)) = spikes(PGHF_times(i):PGHF_times(i+2)-1); 
end
for i = 1:3:length(SGLF_times)
   single_SGLF_spikes((i-1)/3+1, 1:SGLF_times(i+2)-SGLF_times(i)) = spikes(SGLF_times(i):SGLF_times(i+2)-1); 
end
for i = 1:3:length(SGHF_times)
   single_SGHF_spikes((i-1)/3+1, 1:SGHF_times(i+2)-SGHF_times(i)) = spikes(SGHF_times(i):SGHF_times(i+2)-1); 
end

binsize = 1000;
single_PGLF_PSTH_matrix = moving_average_m(single_PGLF_spikes, 3*binsize);
single_PGHF_PSTH_matrix = moving_average_m(single_PGHF_spikes, 3*binsize);
single_SGLF_PSTH_matrix = moving_average_m(single_SGLF_spikes, 3*binsize);
single_SGHF_PSTH_matrix = moving_average_m(single_SGHF_spikes, 3*binsize);

single_PGLF_PSTH = zeros(1, size(single_PGLF_PSTH_matrix,2));
single_PGHF_PSTH = zeros(1, size(single_PGHF_PSTH_matrix,2));
single_SGLF_PSTH = zeros(1, size(single_SGLF_PSTH_matrix,2));
single_SGHF_PSTH = zeros(1, size(single_SGHF_PSTH_matrix,2));

for i = 1:size(single_PGLF_PSTH_matrix,1)
    single_PGLF_PSTH = single_PGLF_PSTH + (1/size(single_PGLF_PSTH_matrix,1))*single_PGLF_PSTH_matrix(i,:);
end
for i = 1:size(PGHF_PSTH_matrix,1)
    single_PGHF_PSTH = single_PGHF_PSTH + (1/size(single_PGHF_PSTH_matrix,1))*single_PGHF_PSTH_matrix(i,:);
end
for i = 1:size(SGLF_PSTH_matrix,1)
    single_SGLF_PSTH = single_SGLF_PSTH + (1/size(single_SGLF_PSTH_matrix,1))*single_SGLF_PSTH_matrix(i,:);
end
for i = 1:size(SGHF_PSTH_matrix,1)
    single_SGHF_PSTH = single_SGHF_PSTH + (1/size(single_SGHF_PSTH_matrix,1))*single_SGHF_PSTH_matrix(i,:);
end

%% compare different events on one neuron

plot(linspace(0, length(single_PGLF_PSTH)/fs, length(single_PGLF_PSTH)), single_PGLF_PSTH*fs); hold on;
plot(linspace(0, length(single_PGHF_PSTH)/fs, length(single_PGHF_PSTH)), single_PGHF_PSTH*fs); hold on;
plot(linspace(0, length(single_SGLF_PSTH)/fs, length(single_SGLF_PSTH)), single_SGLF_PSTH*fs); hold on;
plot(linspace(0, length(single_SGHF_PSTH)/fs, length(single_SGHF_PSTH)), single_SGHF_PSTH*fs); hold on;

        title(['PSTH of event-based spikes (neuron id: ', num2str(neuron_id), ')']); xlabel('task time(s)');
        xline(0.8, 'r-', 'LineWidth', 2); hold on; 
        xline(2.1, 'b-', 'LineWidth', 2); hold on;
        xline(0.4, 'c-', 'LineWidth', 2); hold on;
        xline(SRONSET_avg/fs, 'g-', 'LineWidth', 2);xlim([0, 5.6]); legend({'PGLF PSTH', 'PGHF PSTH', 'SGLF PSTH', 'SGHF PSTH', 'CUE-ON', 'GO-ON', 'WS-ON', 'SR-ON'}, 'Location', 'northeast'); ylabel('neuron firing rate (Hz)');
        
%% ISI histogram
clc;
ISI = [];
for i = 1:size(block.segments{1,1}.spiketrains,2)
   for j = 1:size(block.segments{1,1}.spiketrains{1,i}.times, 2)-1
       ISI = [ISI, block.segments{1,1}.spiketrains{1,i}.times(j+1)-block.segments{1,1}.spiketrains{1,i}.times(j)];
   end
end

%% save ISI vector
fs = 30000;
ISI_vector = ISI / fs;
save('ISI_vector', 'ISI_vector');

%% plot ISI histogram
clc;
load('ISI_vector.mat');
figure;
lambda = 1/mean(ISI_vector);
hist(ISI_vector, 50000); hold on; xlim([0 1.5]); title('ISI histogram'); xlabel('t(s)'); ylabel('population');

%%
mu = mean(ISI_vector);
sigma = var(ISI_vector);
dist = [];
for i=1:size(ISI_vector, 2)
    if(ISI_vector(1,i) < mu + 3 * sigma)
        dist = [dist, ISI_vector(1,i)];
    end
end
histogram(dist)

%%
mean(ISI_vector)
var(ISI_vector)

%%
bin_size = 0.01;
D = zeros(65536);
for i=1:size(ISI_vector, 2)
    j = uint(ISI_vector(i)/bin_size);
    D(j) = D(j)+1;
end
plot(1:65536, D)

%%
figure
histogram(ISI_vector, 100000);
xlim([0 1.5]); title('ISI histogram'); xlabel('t(s)'); ylabel('population');
mu = mean(ISI_vector);
X_exp = 0:0.01:1.5;
pdf = exppdf(X_exp, mu);
hold on
plot(X_exp, pdf, 'g')

%%
figure
histfit(ISI_vector, 50000, 'exponential');
xlim([0 1.5]); title('ISI histogram'); xlabel('t(s)'); ylabel('population');
%pd = fitdist(ISI_vector','gamma');
%mu = mean(ISI_vector);
%X = 0:0.01:1.5;
%pdf = exppdf(X,mu);
%plot(X, pdf*2002091, 'g')

%%
BinEdges_exp = linspace(0,10,2501);
cnt_exp = histcounts(ISI_vector,'BinEdges',BinEdges_exp);
X_exp = 0.002 : 0.004 : 9.998;
bar(X_exp, cnt_exp)
tbl_exp = table(X_exp', cnt_exp');

%%
expfun = @(b,x) b(1) * exp(b(2)*x(:, 1)); 
plot(X_exp, cnt_exp, 'b*', 'LineWidth', 2, 'MarkerSize', 15);

%%
exp0 = [25000, -1];
exp_dist = fitnlm(tbl_exp,expfun, exp0);

%%
bar(X_exp, cnt_exp)
coefficients = exp_dist.Coefficients{:, 'Estimate'}
xlim([0 1.5]); ylim([0 50000]); title('ISI histogram'); xlabel('t(s)'); ylabel('population');
yFitted = coefficients(1) * exp(coefficients(2)*X_exp);
hold on;
plot(X_exp, yFitted, 'r-', 'LineWidth', 2);
grid on;
title('Exponential Distribution Fitting');
legend('ISI histogram', 'Exponential Distribution', 'Location', 'north');

%%

BinEdges_gamma = linspace(0,0.01,1001);
cnt_gamma = histcounts(ISI_vector,'BinEdges',BinEdges_gamma);
X_gamma = 0.000005 : 0.00001 : 0.009995;
tbl_gamma = table(X_gamma', cnt_gamma');
bar(X_gamma*1000, cnt_gamma);
title('Zoomed ISI histogram'); xlabel('t(ms)'); ylabel('population');

%%
gammafun = @(b,x) b(1) .* x(:,1).^(b(2)-1) .* exp(b(3)*x(:, 1));
gamma0 = [121090, 0.00001, -240];
gamma_dist = fitnlm(tbl_gamma, gammafun, gamma0);

%%
bar(X_gamma*1000, cnt_gamma)
coefficients = gamma_dist.Coefficients{:, 'Estimate'};

title('ISI histogram'); xlabel('t(ms)'); ylabel('population');
yFitted = coefficients(1)* X_gamma(:,1).^(coefficients(2)-1) * exp(coefficients(3)*X_gamma);
hold on;
plot(X_gamma*1000, yFitted, 'r-', 'LineWidth', 2);
grid on;
title('Gamma Distribution Fitting');
legend('ISI histogram', 'Gamma Distribution', 'Location', 'north');
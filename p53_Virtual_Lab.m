function p53_Virtual_Lab()
    
    clc; close all; clear;
   
    Damage_Levels = [1, 5, 50]; 
    Labels = {'Threshold (0.1 Gy)', 'Low (0.5 Gy)', 'High (5.0 Gy)'};
    
    N_Trials = 100; 
    
    
    Results = struct('Period', [], 'Amplitude', [], 'NumPeaks', []);
    
    figure('Color','w', 'Name', 'Single Trajectories'); 
    
    
    for i = 1:length(Damage_Levels)
        dose = Damage_Levels(i);
        fprintf('Simulating Condition: %s...\n', Labels{i});
        
        periods_of_responders = [];
        amplitudes_of_responders = [];
        peaks_all_cells = [];
        count_responders = 0;
        
        for n = 1:N_Trials
            [T, X] = run_gillespie(dose);
            
            
            [p_mean, a_mean, n_peaks, T_smooth, p53_smooth] = analyze_trace(T, X);
            
            peaks_all_cells = [peaks_all_cells, n_peaks];
            
            
            if n_peaks >= 3
                count_responders = count_responders + 1;
                periods_of_responders = [periods_of_responders, p_mean];
                amplitudes_of_responders = [amplitudes_of_responders, a_mean];
            end
            
            
            if n == 1
                subplot(3,1,i);
                plot(T/3600, X(:,1)+X(:,3), 'Color', [0.8 0.8 0.8]); hold on;
                plot(T_smooth/3600, p53_smooth, 'r', 'LineWidth', 2);
                yline(160, 'b--'); % Show the threshold line
                title([Labels{i}, ' (Single Trial)']);
                ylabel('p53 Level');
                xlim([0 48]); ylim([0 600]);
            end
        end     
        
        
        Results(i).FractionOscillating = (count_responders / N_Trials) * 100; 
        
        if count_responders > 0
            Results(i).Period = mean(periods_of_responders);
            Results(i).Amplitude = mean(amplitudes_of_responders);
            
            Results(i).AvgPeaks = mean(peaks_all_cells(peaks_all_cells >= 3)); 
        else
            Results(i).Period = 0;
            Results(i).Amplitude = 0;
            Results(i).AvgPeaks = 0;
        end
    end
    
    
    figure('Color','w', 'Name', 'Digital Signaling Analysis');
    
    % FRACTION OF RESPONDERS (Probability)
    subplot(1,3,1);
    bar([Results.FractionOscillating], 'FaceColor', [0.2 0.6 0.8]);
    xticklabels({'0.1 Gy', '0.5 Gy', '5.0 Gy'});
    ylabel('% of Oscillating Cells');
    title('Probability (Responders)');
    ylim([0 100]);
    grid on;
    
    % PEAKS
    subplot(1,3,2);
    bar([Results.AvgPeaks], 'FaceColor', [0.8 0.4 0.2]);
    xticklabels({'0.1 Gy', '0.5 Gy', '5.0 Gy'});
    ylabel('Avg Peaks (Responders Only)');
    title('Duration (Pulse Count)');
    grid on;
    
    % PERIOD
    subplot(1,3,3);
    bar([Results.Period], 'FaceColor', [0.4 0.8 0.4]);
    xticklabels({'0.1 Gy', '0.5 Gy', '5.0 Gy'});
    ylabel('Avg Period (Hours)');
    title('Robustness (Clock Speed)');
    ylim([0 8]);
    grid on;

    fprintf('\n--- FINAL RESULTS ---\n');
    for i = 1:3
        fprintf('%s -> Responders: %.0f%% | Peaks (Responders): %.1f | Period: %.2f h | Amp: %.0f\n', ...
            Labels{i}, Results(i).FractionOscillating, Results(i).AvgPeaks, Results(i).Period, Results(i).Amplitude);
    end
end


function [avg_period, avg_amp, num_peaks, T_reg, p53_smooth] = analyze_trace(T, X)
    Total_p53 = X(:,1) + X(:,3);
    
    % Interpolate
    T_reg = linspace(min(T), max(T), 5000)';
    p53_reg = interp1(T, Total_p53, T_reg, 'linear');
    window_size = 100;
    kernel = ones(window_size, 1) / window_size;
    p53_smooth = conv(p53_reg, kernel, 'same');
    
 
    pks = []; locs = [];
    

    min_height = 160; 
    
    min_dist_indices = 300; 
    last_loc = -min_dist_indices;
    
    for i = 2:length(p53_smooth)-1
        if p53_smooth(i) > p53_smooth(i-1) && p53_smooth(i) > p53_smooth(i+1)
            if p53_smooth(i) > min_height 
                if (i - last_loc) > min_dist_indices
                    pks = [pks, p53_smooth(i)];
                    locs = [locs, T_reg(i)];
                    last_loc = i;
                else
                    if p53_smooth(i) > pks(end)
                        pks(end) = p53_smooth(i);
                        locs(end) = T_reg(i);
                        last_loc = i;
                    end
                end
            end
        end
    end
                        
    num_peaks = length(pks);
    
    if num_peaks > 1
        diffs = diff(locs) / 3600; 
        avg_period = mean(diffs);
        avg_amp = mean(pks);
    else
        avg_period = NaN; 
        if num_peaks == 1, avg_amp = pks(1); else, avg_amp = NaN; end
    end
end

%GILLESPIE SIMULATION 
function [T_hist, X_hist] = run_gillespie(IR_dose)

    ksynMdm2 = 4.95e-4; kdegMdm2 = 4.33e-4;
    ksynp53 = 0.078;    kdegp53 = 8.25e-4;
    kbinMdm2p53 = 11.55e-4; krelMdm2p53 = 1.155e-5;
    ksynMdm2mRNA = 1.0e-4; kdegMdm2mRNA = 1.0e-4;
    kbinARFMdm2 = 1.0e-2; kdegARFMdm2 = 1.0e-3;
    
    
    kdegARF = 2.0e-4;    
    kactARF = 3.3e-5;
    kdam = 0.8;          
    krepair = 1.0e-3;    

    tspan = 48*3600;     
    
    % Initial State: [p53, Mdm2, Mdm2_p53, Mdm2_mRNA, ARF, ARF_Mdm2, damDNA]
    x = [5, 5, 95, 0, 0, 0, 0]; 
    t = 0;
    
    ir_time = 3600;
    

    T_hist = zeros(100000, 1);
    X_hist = zeros(100000, 7);
    count = 1; T_hist(1) = 0; X_hist(1,:) = x;
    
    while t < tspan
        % Calc Propensities
        if t >= ir_time && t < ir_time+60, cur_IR = IR_dose; else, cur_IR = 0; end
        
        a = zeros(14,1);
        a(1) = ksynMdm2 * x(4); 
        a(2) = ksynMdm2mRNA * x(1);
        a(3) = kdegMdm2mRNA * x(4);
        a(4) = kdegMdm2 * x(2);
        a(5) = ksynp53;
        a(6) = kdegp53 * x(3);
        a(7) = kbinMdm2p53 * x(1) * x(2);
        a(8) = krelMdm2p53 * x(3);
        a(9) = kdam * cur_IR;
        a(10) = krepair * x(7);
        a(11) = kactARF * x(7);
        a(12) = kbinARFMdm2 * x(5) * x(2);
        a(13) = kdegARFMdm2 * x(6);
        a(14) = kdegARF * x(5);
        
        a0 = sum(a);
        if a0 == 0, break; end
        
        tau = (1/a0)*log(1/rand);
        
        r2 = rand * a0; sum_a = 0; mu = 0;
        for k = 1:14
            sum_a = sum_a + a(k);
            if sum_a >= r2, mu = k; break; end
        end
        
        if mu==1, x(2)=x(2)+1;
        elseif mu==2, x(4)=x(4)+1;
        elseif mu==3, x(4)=x(4)-1;
        elseif mu==4, x(2)=x(2)-1;
        elseif mu==5, x(1)=x(1)+1;
        elseif mu==6, x(3)=x(3)-1; x(2)=x(2)+1;
        elseif mu==7, x(1)=x(1)-1; x(2)=x(2)-1; x(3)=x(3)+1;
        elseif mu==8, x(3)=x(3)-1; x(1)=x(1)+1; x(2)=x(2)+1;
        elseif mu==9, x(7)=x(7)+1;
        elseif mu==10, x(7)=x(7)-1;
        elseif mu==11, x(5)=x(5)+1;
        elseif mu==12, x(5)=x(5)-1; x(2)=x(2)-1; x(6)=x(6)+1;
        elseif mu==13, x(6)=x(6)-1; x(5)=x(5)+1;
        elseif mu==14, x(5)=x(5)-1;
        end
        
        t = t + tau;
        count = count + 1;
        if count > length(T_hist) 
            T_hist = [T_hist; zeros(50000,1)];
            X_hist = [X_hist; zeros(50000,7)];
        end
        T_hist(count) = t; X_hist(count,:) = x;
    end
    T_hist = T_hist(1:count); X_hist = X_hist(1:count,:);
end
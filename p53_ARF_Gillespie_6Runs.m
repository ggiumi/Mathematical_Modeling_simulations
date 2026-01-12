function p53_ARF_Gillespie_6Runs()
   
    
    clc; close all; clear;

    ksynMdm2      = 4.95e-4;
    kdegMdm2      = 4.33e-4;   
    ksynp53       = 0.078;
    kdegp53       = 8.25e-4;
    kbinMdm2p53   = 11.55e-4;
    krelMdm2p53   = 1.155e-5;
    ksynMdm2mRNA  = 1.0e-4;
    kdegMdm2mRNA  = 1.0e-4;
    kbinARFMdm2   = 1.0e-2;
    kdegARFMdm2   = 1.0e-3;
    kdegARF       = 1.0e-4;
    kactARF       = 3.3e-5;
    kdam          = 0.08;
    krepair       = 2.0e-5;
    
    tspan = [0, 30*3600]; % 30 hours
   
    IR_dose = 25;      
    irradiation_time = 3600; 
    
    figure('Color','w', 'Name', 'p53-Mdm2 Stochastic Variability (6 Runs)');
    
    num_simulations = 6;
    
    for sim = 1:num_simulations
        fprintf('Running Simulation %d of %d...\n', sim, num_simulations);
        

        x = [5, 5, 95, 0, 0, 0, 0]; 
        t = 0;
        
        max_steps = 1000000;
        T_hist = zeros(max_steps, 1);
        X_hist = zeros(max_steps, 7);
        
        T_hist(1) = t;
        X_hist(1,:) = x;
        count = 1;
        
        while t < tspan(2)
            
            p53 = x(1); Mdm2 = x(2); Mdm2_p53 = x(3); Mdm2_mRNA = x(4);
            ARF = x(5); ARF_Mdm2 = x(6); damDNA = x(7);
            
          
            if t >= irradiation_time && t < irradiation_time + 60
                 Current_IR = IR_dose;
            else
                 Current_IR = 0;
            end
            
            a = zeros(14, 1);
            a(1)  = ksynMdm2 * Mdm2_mRNA;      
            a(2)  = ksynMdm2mRNA * p53;        
            a(3)  = kdegMdm2mRNA * Mdm2_mRNA;  
            a(4)  = kdegMdm2 * Mdm2;           
            a(5)  = ksynp53;                   
            a(6)  = kdegp53 * Mdm2_p53;        
            a(7)  = kbinMdm2p53 * p53 * Mdm2;  
            a(8)  = krelMdm2p53 * Mdm2_p53;    
            a(9)  = kdam * Current_IR;         
            a(10) = krepair * damDNA;          
            a(11) = kactARF * damDNA;          
            a(12) = kbinARFMdm2 * ARF * Mdm2;  
            a(13) = kdegARFMdm2 * ARF_Mdm2;    
            a(14) = kdegARF * ARF;             
            
            a0 = sum(a);
            if a0 == 0, break; end
            
            tau = (1/a0) * log(1/rand);
            r2 = rand * a0;
            sum_a = 0; mu = 0;
            for k=1:14
                sum_a = sum_a + a(k);
                if sum_a >= r2, mu = k; break; end
            end
            

            if mu == 1, x(2)=x(2)+1;
            elseif mu == 2, x(4)=x(4)+1;
            elseif mu == 3, x(4)=x(4)-1;
            elseif mu == 4, x(2)=x(2)-1;
            elseif mu == 5, x(1)=x(1)+1;
            elseif mu == 6, x(3)=x(3)-1; x(2)=x(2)+1;
            elseif mu == 7, x(1)=x(1)-1; x(2)=x(2)-1; x(3)=x(3)+1;
            elseif mu == 8, x(3)=x(3)-1; x(1)=x(1)+1; x(2)=x(2)+1;
            elseif mu == 9, x(7)=x(7)+1;
            elseif mu == 10, x(7)=x(7)-1;
            elseif mu == 11, x(5)=x(5)+1;
            elseif mu == 12, x(5)=x(5)-1; x(2)=x(2)-1; x(6)=x(6)+1;
            elseif mu == 13, x(6)=x(6)-1; x(5)=x(5)+1;
            elseif mu == 14, x(5)=x(5)-1;
            end
            
            t = t + tau;
            count = count + 1;
            
            if count > size(T_hist,1)
                 T_hist = [T_hist; zeros(100000,1)];
                 X_hist = [X_hist; zeros(100000,7)];
            end
            T_hist(count) = t;
            X_hist(count,:) = x;
        end
        
        T_hist = T_hist(1:count);
        X_hist = X_hist(1:count, :);
        
        subplot(2, 3, sim); 
        
        Total_p53 = X_hist(:,1) + X_hist(:,3); 
        Total_Mdm2 = X_hist(:,2) + X_hist(:,3) + X_hist(:,6);
        
        plot(T_hist/3600, Total_p53, 'b', 'LineWidth', 1); hold on;
        plot(T_hist/3600, Total_Mdm2, 'r', 'LineWidth', 1);
        
        xlim([0 30]); ylim([0 600]); 
  
        if sim > 3, xlabel('Time (h)'); end
        if mod(sim,3)==1, ylabel('Molecules'); end
        grid on;
    end
    sgtitle('Stochastic Simulation of p53-Mdm2 (ARF Model)', 'FontSize', 14, 'FontWeight', 'bold');

    
    subplot(2,3,1); legend('Total p53', 'Total Mdm2');
    
end
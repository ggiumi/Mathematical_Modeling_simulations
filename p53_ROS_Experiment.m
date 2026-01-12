function p53_ROS_experiment()

    
    clc; close all; clear;


    P.ksynMdm2 = 4.95e-4; 
    P.kdegMdm2 = 4.33e-4 ;      
    P.ksynp53 = 0.078;    
    P.kdegp53 = 8.25e-4 ;       
    P.kbinMdm2p53 = 11.55e-4; 
    P.krelMdm2p53 = 1.155e-5;
    P.ksynMdm2mRNA = 1.0e-4; 
    P.kdegMdm2mRNA = 1.0e-4;
    P.kbinARFMdm2 = 1.0e-2; 
    P.kdegARFMdm2 = 1.0e-3 ;    
    P.kdegARF = 1.0e-4*5 ;         
    P.kactARF = 3.3e-5;
    
   
    P.ksynROS = 0.001;     
    P.kdamROS = 3.0e-6;    
    
    P.krepair_LOW = 1.0e-5;    
    P.krepair_HIGH = 1.0e-1;   
    P.repair_time = 80 * 3600; 

    tspan = [0, 120*3600]; % 120 Hours
    
   
    figure('Color','w', 'Name', 'Figure 8 Corrected', 'Position', [50 50 1200 800]);

    num_sims = 6;
    
    for i = 1:num_sims
        fprintf('Running Simulation %d of %d...\n', i, num_sims);
        
        % Initial Conditions [p53, Mdm2, Complex, mRNA, ARF, ARF_Mdm2, damDNA, ROS]
        x0 = [10, 10, 100, 0, 0, 0, 0, 0]; 
        
        
        [T, X] = run_gillespie_ROS(P, tspan, x0);
        
        
        subplot(2, 3, i);
        
        
        Total_p53 = X(:,1) + X(:,3);        
        Total_Mdm2 = X(:,2) + X(:,3) + X(:,6); 
        damDNA = X(:,7);
        
        
        plot(T/3600, Total_p53, 'g', 'LineWidth', 1); hold on;
        plot(T/3600, Total_Mdm2, 'r', 'LineWidth', 1);
        plot(T/3600, damDNA, 'k', 'LineWidth', 1.5); % Black line
        
        
        xlim([0 120]);
        ylim([0 400]); 
        xlabel('Time (hours)');
        ylabel('Number of proteins');
        
        
        if i == 1
            legend('p53', 'Mdm2', 'damaged DNA', 'Location', 'NorthWest');
        end
    end
end

function [T_hist, X_hist] = run_gillespie_ROS(P, tspan, x0)
    
    ksynMdm2 = P.ksynMdm2; kdegMdm2 = P.kdegMdm2;
    ksynp53 = P.ksynp53; kdegp53 = P.kdegp53;
    kbinMdm2p53 = P.kbinMdm2p53; krelMdm2p53 = P.krelMdm2p53;
    ksynMdm2mRNA = P.ksynMdm2mRNA; kdegMdm2mRNA = P.kdegMdm2mRNA;
    kbinARFMdm2 = P.kbinARFMdm2; kdegARFMdm2 = P.kdegARFMdm2;
    kdegARF = P.kdegARF; kactARF = P.kactARF;
    
    ksynROS = P.ksynROS; kdamROS = P.kdamROS;
    krepair_LOW = P.krepair_LOW; krepair_HIGH = P.krepair_HIGH;
    repair_time = P.repair_time;
    
    x = x0; t = 0;
    

    chunk = 500000;
    T_hist = zeros(chunk, 1);
    X_hist = zeros(chunk, 8);
    count = 1; T_hist(1)=0; X_hist(1,:)=x;
    
    while t < tspan(2)
        if t < repair_time
            curr_krepair = krepair_LOW;
        else
            curr_krepair = krepair_HIGH;
        end
        
        % Propensities
        a = zeros(15,1);
        a(1) = ksynMdm2 * x(4); 
        a(2) = ksynMdm2mRNA * x(1);
        a(3) = kdegMdm2mRNA * x(4);
        a(4) = kdegMdm2 * x(2);
        a(5) = ksynp53;
        a(6) = kdegp53 * x(3);
        a(7) = kbinMdm2p53 * x(1) * x(2);
        a(8) = krelMdm2p53 * x(3);
        a(9) = kdamROS * x(8);         
        a(10)= curr_krepair * x(7);    
        a(11)= kactARF * x(7);
        a(12)= kbinARFMdm2 * x(5) * x(2);
        a(13)= kdegARFMdm2 * x(6);
        a(14)= kdegARF * x(5);
        a(15)= ksynROS;                
        
        a0 = sum(a); if a0 == 0, break; end
        
        tau = (1/a0)*log(1/rand);
        r2 = rand * a0; sum_a = 0; mu = 0;
        for k=1:15, sum_a=sum_a+a(k); if sum_a>=r2, mu=k; break; end; end
        
     
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
        elseif mu==15, x(8)=x(8)+1;
        end
        
        t = t + tau; count = count + 1;
        
        if count > size(T_hist,1)
            T_hist = [T_hist; zeros(chunk,1)];
            X_hist = [X_hist; zeros(chunk,8)];
        end
        T_hist(count) = t; X_hist(count,:) = x;
    end
    T_hist = T_hist(1:count); X_hist = X_hist(1:count,:);
    sgtitle('Stochastic Simulation of p53-Mdm2 with ROS-induced damage', 'FontSize', 14, 'FontWeight', 'bold');
end
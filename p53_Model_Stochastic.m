function p53_Model_integrated
 
    clear; clc; close all;
    
    %% --- 1. CONFIGURATION ---
    INCLUDE_WIP1   = true;  
    INCLUDE_MIR192 = true;  
    
    t_max_hours = 40;       
    t_max = t_max_hours * 3600;
    
    %% --- 2. INITIALIZATION ---
    x = zeros(14, 1);
    x(1) = 10; x(2) = 10; x(6) = 200; x(10)= 10;
    
    max_steps = 2e6; 
    T_hist = zeros(max_steps, 1);
    X_hist = zeros(max_steps, 14);
    count = 1; X_hist(1,:) = x; t = 0;
    
    %% --- 3. PARAMETERS ---
    k_syn_p53_mRNA = 0.001; k_deg_p53_mRNA = 0.0001;
    k_syn_p53 = 0.006; k_deg_p53 = 0.000825;
    k_bin = 1.155e-3; k_rel = 1.155e-5;
    k_syn_mdm2_rna = 0.0001; k_deg_mdm2_rna = 0.0001;
    k_syn_mdm2 = 0.000495; k_deg_mdm2 = 0.000433;
    k_deg_mdm2_atm = 0.0004;
    k_dam = 0.08; k_repair = 2e-5;
    k_act_atm = 1e-4; k_inact_atm = 5e-4;
    k_phos_p53 = 5e-4; k_dephos_p53 = 0.5;
    k_phos_mdm2 = 2.0; k_dephos_mdm2= 0.5;
    k_syn_wip1_rna = 0.0001; k_deg_wip1_rna = 0.0001;
    k_syn_wip1 = 0.0005; k_deg_wip1 = 0.0005;
    k_wip1_act = 0.01;
    k_syn_mir_pre = 0.0001; k_deg_mir_pre = 0.0001;
    k_mat_mir = 0.002; k_deg_mir_mat = 0.0001;
    
    % Fixed Parameter
    k_mir_action = 0.00001; 
    
    %% --- 4. SIMULATION ---
    fprintf('Running Simulation...\n');
    while t < t_max
        p53=x(1); Mdm2=x(2); Mdm2_P=x(3); Mdm2_R=x(4); Cmplx=x(5);
        ATMI=x(6); ATMA=x(7); Dam=x(8); p53_P=x(9); p53_R=x(10);
        Wip1=x(12); miR_Mat=x(14); Wip1_R=x(11); miR_Pre=x(13);
        
        IR = 0; if t > 3600 && t < 3660, IR = 25; end
        
        a = zeros(30, 1);
        a(1) = k_syn_p53_mRNA;
        a(2) = k_deg_p53_mRNA * p53_R;
        a(3) = k_syn_p53 * p53_R;
        a(4) = k_bin * p53 * Mdm2;
        a(5) = k_rel * Cmplx;
        a(6) = k_deg_p53 * Cmplx;
        a(7) = k_dam * IR;
        a(8) = k_repair * Dam;
        a(9) = k_act_atm * Dam * ATMI;
        
        r_inact = k_inact_atm; if INCLUDE_WIP1, r_inact = r_inact + k_wip1_act * Wip1; end
        a(10) = r_inact * ATMA;
        
        a(11) = k_phos_p53 * p53 * ATMA;
        a(12) = k_phos_mdm2 * Mdm2 * ATMA;
        
        r_dp53 = k_dephos_p53; if INCLUDE_WIP1, r_dp53 = r_dp53 + k_wip1_act * Wip1; end
        a(13) = r_dp53 * p53_P;
        a(14) = k_dephos_mdm2 * Mdm2_P;
        
        a(15) = k_syn_mdm2_rna * p53_P; 
        rate_mdm2_deg = k_deg_mdm2_rna;
        if INCLUDE_MIR192, rate_mdm2_deg = rate_mdm2_deg + (k_mir_action * miR_Mat); end
        a(16) = rate_mdm2_deg * Mdm2_R;
        
        a(17) = k_syn_mdm2 * Mdm2_R;
        a(18) = k_deg_mdm2 * Mdm2;
        a(19) = k_deg_mdm2_atm * Mdm2_P;
        a(20) = k_deg_p53 * p53_P;
        
        if INCLUDE_WIP1
            a(21) = k_syn_wip1_rna * p53_P; a(22) = k_deg_wip1_rna * Wip1_R;
            a(23) = k_syn_wip1 * Wip1_R; a(24) = k_deg_wip1 * Wip1;
        end
        if INCLUDE_MIR192
            a(25) = k_syn_mir_pre * p53_P; a(26) = k_deg_mir_pre * miR_Pre;
            a(27) = k_mat_mir * miR_Pre;   a(28) = k_deg_mir_mat * miR_Mat;
        end
        
        a0 = sum(a); if a0 == 0, break; end
        t = t + (-log(rand) / a0);
        r = rand * a0;
        idx = find(cumsum(a) >= r, 1);
        
        % State Update
        switch idx
            case 1, x(10)=x(10)+1; case 2, x(10)=x(10)-1; case 3, x(1)=x(1)+1;
            case 4, x(1)=x(1)-1; x(2)=x(2)-1; x(5)=x(5)+1;
            case 5, x(1)=x(1)+1; x(2)=x(2)+1; x(5)=x(5)-1;
            case 6, x(5)=x(5)-1; x(2)=x(2)+1;
            case 7, x(8)=x(8)+1; case 8, x(8)=x(8)-1;
            case 9, x(6)=x(6)-1; x(7)=x(7)+1; case 10,x(6)=x(6)+1; x(7)=x(7)-1;
            case 11,x(1)=x(1)-1; x(9)=x(9)+1; case 12,x(2)=x(2)-1; x(3)=x(3)+1;
            case 13,x(9)=x(9)-1; x(1)=x(1)+1; case 14,x(3)=x(3)-1; x(2)=x(2)+1;
            case 15,x(4)=x(4)+1; case 16,x(4)=x(4)-1; 
            case 17,x(2)=x(2)+1; case 18,x(2)=x(2)-1; case 19,x(3)=x(3)-1;
            case 20,x(9)=x(9)-1;
            case 21,x(11)=x(11)+1; case 22,x(11)=x(11)-1;
            case 23,x(12)=x(12)+1; case 24,x(12)=x(12)-1;
            case 25,x(13)=x(13)+1; case 26,x(13)=x(13)-1; 
            case 27,x(13)=x(13)-1; x(14)=x(14)+1; case 28,x(14)=x(14)-1; 
        end
        if count < max_steps, count=count+1; T_hist(count)=t; X_hist(count,:)=x; end
    end
    T_hist = T_hist(1:count); X_hist = X_hist(1:count,:); T_hrs = T_hist/3600;
    
    %% --- 5. PLOTTING (SINGLE AXIS) ---
    Total_p53  = X_hist(:,1) + X_hist(:,9) + X_hist(:,5);
    Total_Mdm2 = X_hist(:,2) + X_hist(:,3) + X_hist(:,5);
    
    figure('Color','w', 'Name', 'Single_Axis_Result');
    
    subplot(4,1,1);
    % Standard Plot: Both on same axis
    plot(T_hrs, Total_p53, 'b', 'LineWidth', 1.5); hold on;
    plot(T_hrs, Total_Mdm2, 'r', 'LineWidth', 1.5);
    
    title('Total p53 (Blue) vs Total Mdm2 (Red)');
    ylabel('Molecules');
    legend('Total p53', 'Total Mdm2');
    grid on;
    
    subplot(4,1,2); plot(T_hrs, X_hist(:,7), 'c', 'LineWidth', 1); title('ATM Active'); grid on;
    subplot(4,1,3); plot(T_hrs, X_hist(:,12), 'm', 'LineWidth', 1); title('Wip1 Protein'); grid on;
    subplot(4,1,4); plot(T_hrs, X_hist(:,14), 'Color', [0 0.6 0], 'LineWidth', 1); title('miR-192 Mature'); xlabel('Time (Hours)'); grid on;
end
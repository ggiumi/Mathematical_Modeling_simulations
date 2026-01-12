function p53_Integrated_SanityCheck
    
    clear; clc; close all;
    
    %% --- 1. EXPERIMENTAL TOGGLES ---
    % 1. Start with FALSE/FALSE to see the clean Proctor Oscillator.
    % 2. Turn on WIP1 -> Oscillations should die out faster (Termination).
    % 3. Turn on MIR192 -> Oscillations should change shape (Robustness).
    
    ENABLE_WIP1   = true; 
    ENABLE_MIR192 = true;
    
    % Simulation Settings
    t_max_hours = 24;       
    t_max = t_max_hours * 3600;
    
    %% --- 2. INITIALIZATION (Steady State) ---
    % We start with non-zero values to prevent immediate system crash.
    x = zeros(14, 1);
    
    % Core Species
    x(1) = 20;   % p53
    x(2) = 10;   % Mdm2
    x(5) = 100;  % p53-Mdm2 Complex (The reservoir)
    x(6) = 100;  % ATM (Inactive)
    
    % RNA
    x(4) = 10;   % Mdm2 mRNA
    x(10)= 10;   % p53 mRNA
    
    %% --- 3. PARAMETERS (Calibrated for Compatibility) ---
    
    % --- CORE (Proctor-like) ---
    k_prod_p53      = 0.05;    % Protein synthesis
    k_deg_p53       = 0.001;   % Basal degradation
    
    k_bind          = 0.01;    % p53 + Mdm2 binding
    k_rel           = 0.0001;  % Dissociation
    k_deg_cmplx     = 0.005;   % Degradation of p53 in complex
    
    % Mdm2 Production (The Critical Link)
    k_prod_mdm2_rna = 0.001;   % Basal
    k_ind_mdm2_rna  = 0.01;    % Induced by p53_P
    k_deg_mdm2_rna  = 0.001;
    k_prod_mdm2     = 0.05;
    k_deg_mdm2      = 0.002;
    
    % ATM & Damage (The Driver)
    k_dam_on        = 0.5;
    k_repair        = 0.0001;  % Slow repair = Sustained oscillations
    k_atm_act       = 0.002;
    k_atm_inact     = 0.001;
    
    % Phosphorylation
    k_phos_p53      = 0.005;
    k_dephos_p53    = 0.05;
    
    % --- THE LOOPS (Relative Strengths) ---
    
    % Wip1 (The Terminator)
    % Produced by p53_P. Destroys p53_P and Inactivates ATM.
    k_prod_wip1     = 0.005;
    k_deg_wip1      = 0.001;
    wip1_strength   = 0.005;   % How hard does Wip1 hit?
    
    % miR-192 (The Stabilizer)
    % Produced by p53_P. Destroys Mdm2 mRNA.
    k_prod_mir      = 0.002;
    k_deg_mir       = 0.0005;
    mir_strength    = 0.002;   % How hard does miR hit Mdm2?
    
    %% --- 4. SIMULATION ---
    t = 0;
    count = 1;
    max_steps = 1e6;
    T_hist = zeros(max_steps, 1);
    X_hist = zeros(max_steps, 14);
    X_hist(1,:) = x;
    
    fprintf('Running integrated simulation...\n');
    
    while t < t_max
        
        % Unpack
        p53=x(1); Mdm2=x(2); Mdm2_R=x(4); Cmplx=x(5);
        ATMI=x(6); ATMA=x(7); Dam=x(8); p53_P=x(9); p53_R=x(10);
        Wip1=x(12); miR=x(14);
        
        % Damage Event: Irradiation at t=1h
        IR = 0;
        if t > 3600 && t < 3700, IR = 20; end
        
        % Propensities
        a = zeros(20,1);
        
        % 1. p53 Production/Degradation
        a(1) = 0.005;               % mRNA const
        a(2) = k_prod_p53 * p53_R;
        a(3) = k_deg_p53 * p53;
        
        % 2. Interaction
        a(4) = k_bind * p53 * Mdm2;
        a(5) = k_rel * Cmplx;
        a(6) = k_deg_cmplx * Cmplx;
        
        % 3. ATM / Damage
        a(7) = k_dam_on * IR;
        a(8) = k_repair * Dam;
        a(9) = k_atm_act * Dam * ATMI;
        
        % Wip1 Loop: Increases ATM inactivation
        r_atm_off = k_atm_inact;
        if ENABLE_WIP1, r_atm_off = r_atm_off + (wip1_strength * Wip1); end
        a(10)= r_atm_off * ATMA;
        
        % 4. Phosphorylation
        a(11) = k_phos_p53 * p53 * ATMA;
        
        % Wip1 Loop: Increases p53 dephosphorylation
        r_dp53 = k_dephos_p53;
        if ENABLE_WIP1, r_dp53 = r_dp53 + (wip1_strength * Wip1); end
        a(12) = r_dp53 * p53_P;
        
        % 5. Mdm2 Cycle
        a(13) = k_prod_mdm2_rna + (k_ind_mdm2_rna * p53_P); % Induced by p53_P
        
        % miR-192 Loop: Increases Mdm2 mRNA degradation
        r_mdm2_deg = k_deg_mdm2_rna;
        if ENABLE_MIR192, r_mdm2_deg = r_mdm2_deg + (mir_strength * miR); end
        a(14) = r_mdm2_deg * Mdm2_R;
        
        a(15) = k_prod_mdm2 * Mdm2_R;
        a(16) = k_deg_mdm2 * Mdm2;
        
        % 6. Wip1 & miR Dynamics
        if ENABLE_WIP1
            a(17) = k_prod_wip1 * p53_P;
            a(18) = k_deg_wip1 * Wip1;
        end
        
        if ENABLE_MIR192
            a(19) = k_prod_mir * p53_P;
            a(20) = k_deg_mir * miR;
        end
        
        % Gillespie Step
        a0 = sum(a);
        if a0 == 0, break; end
        
        t = t + -log(rand)/a0;
        r = rand * a0;
        idx = find(cumsum(a) >= r, 1);
        
        switch idx
            case 1, x(10)=x(10); % Dummy mRNA maintenance
            case 2, x(1)=x(1)+1;
            case 3, x(1)=x(1)-1;
            case 4, x(1)=x(1)-1; x(2)=x(2)-1; x(5)=x(5)+1;
            case 5, x(1)=x(1)+1; x(2)=x(2)+1; x(5)=x(5)-1;
            case 6, x(5)=x(5)-1; x(2)=x(2)+1; % Deg p53, Recycle Mdm2
            case 7, x(8)=x(8)+1;
            case 8, x(8)=x(8)-1;
            case 9, x(6)=x(6)-1; x(7)=x(7)+1;
            case 10,x(6)=x(6)+1; x(7)=x(7)-1;
            case 11,x(1)=x(1)-1; x(9)=x(9)+1;
            case 12,x(9)=x(9)-1; x(1)=x(1)+1;
            case 13,x(4)=x(4)+1;
            case 14,x(4)=x(4)-1;
            case 15,x(2)=x(2)+1;
            case 16,x(2)=x(2)-1;
            case 17,x(12)=x(12)+1;
            case 18,x(12)=x(12)-1;
            case 19,x(14)=x(14)+1;
            case 20,x(14)=x(14)-1;
        end
        
        if count < max_steps
            count = count + 1;
            T_hist(count) = t;
            X_hist(count,:) = x;
        end
    end
    
    % Plotting
    T_hist = T_hist(1:count)/3600;
    X_hist = X_hist(1:count,:);
    
    figure('Color','w'); 
    
    subplot(3,1,1);
    yyaxis left; plot(T_hist, X_hist(:,1)+X_hist(:,9)+X_hist(:,5), 'b'); ylabel('Total p53');
    yyaxis right; plot(T_hist, X_hist(:,2)+X_hist(:,5), 'r'); ylabel('Total Mdm2');
    title('p53 (Blue) vs Mdm2 (Red)');
    
    subplot(3,1,2);
    plot(T_hist, X_hist(:,7), 'c', 'LineWidth', 1.5); title('Active ATM (Damage)');
    
    subplot(3,1,3);
    plot(T_hist, X_hist(:,12), 'm', T_hist, X_hist(:,14), 'g', 'LineWidth', 1.5);
    legend('Wip1','miR-192'); title('Regulators');
end
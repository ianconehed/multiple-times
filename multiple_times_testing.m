clear all
close all

num_patterns = 3;
HD_q = zeros(1,ceil(50/num_patterns)*num_patterns);
for q = 1:1%ceil(50/num_patterns)
%     rng(9)
    rng;

    %Number of units
    N = 50; %number of excitatory units
    n = N; %number of inhibitory units

    t_reward2 = [300, 400, 500, 300, 500, 250, 550, 200, 600, 275, 525, 250, 550];
%     t_reward2 = (4/4)*[400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400];
    t_reward = t_reward2(1:num_patterns);


    %Time parameters
    dt = 1; %time step
    t_total = max(t_reward)+251; %time of trial
    time_t = 1:dt:t_total; %time step array
    tau_exp = 20;
    %Time constants
    tau_u = 50; %excitatory time constant
    tau_v = 10; %inhibitory time constant
    tau_t = 5000;
    tau_a = 5000;
    %Activations and firing rate initialization

    %Activation function parameters
    u_c = 2; %critical synpatic activation
    v = 2;
    theta_phi = 0;

    %Learning parameters
    eta_rec = .00006; %learning rate
    eta_inh = .002;
    eta_a = 40;
    A_max = .5;
    eta_t = 10;
    T_max = 3;
    %input parameters

    I_strength = 5;
    I_ext = zeros(N,t_total);
    I_ext(:,1:100) = I_strength;


    f = 0.2;
    chi = zeros(N,num_patterns);
    for i = 1:num_patterns
        temp = randperm(N,f*N);
        chi(temp,i) = 1;
    end
    patterns = chi;
    crosstalk = zeros(N,num_patterns);
    
    crosstalk(1,1) = (1/N)*(patterns(1,:)*patterns'*patterns(:,1) - patterns(1,1)*patterns(:,1)'*patterns(:,1));
    crosstalk(2,2) = (1/N)*(patterns(1,:)*patterns'*patterns(:,1)) - patterns(1,1);
    crosstalk = (1/N)*(patterns*(patterns')*patterns) - patterns;

    
    
    
    chi = repmat(chi,[1 1 t_total]);
    % chi(1:N/3,2,:) = 0;
    % chi(2*N/3+1:N,1,:) = 0;
    % chi(N/3:N,:,:) = 0;
    A_t = chi;
    chi(:,:,101:t_total) = 0;
    chi = chi*I_strength;

    chi2 = zeros(N,num_patterns);
    for i = 1:num_patterns
        temp = randperm(N,f*N);
        chi2(temp,i) = 1;
    end
    chi2 = repmat(chi2,[1 1 t_total]);
    % chi(1:N/3,2,:) = 0;
    % chi(2*N/3+1:N,1,:) = 0;
    % chi2(N/3:N,:,:) = 0;
    chi2(:,:,101:t_total) = 0;
    chi2 = chi2*I_strength;


    % if learning == 0
    %     I_ext(:,1:100) = I_strength;
    % end

    W_ij = (2*randn(N,N) + .1)/N;
    W_ij = W_ij.*(W_ij>0);
%     W_ij = W_ij - diag(diag(W_ij));
    W_ij = 0*W_ij;

    P_ik = .5*diag(ones(1,N));

    M_ki = .1*rand(N,N)/(N);

    % A_t2 = ones(N,num_patterns,t_total);
    A_t2 = .5*A_t;
    amp = 6.5;


    c_max = 200;
    lambda = .2;
    thresh = 3;

    thresh2 = 0*eta_rec;
    %

    thresh3 = 0.2;
    %% Main Program
    [u_it, v_kt,rate_uit,rate_vkt,T, A_j] = deal(zeros(N,t_total));
    [del_ui] = deal(zeros(N,1));
    rate_uit_ex = cell(1,num_patterns);
    rate_uit_af = zeros(N,t_total,num_patterns);
    ham_ful = zeros(num_patterns,c_max);

    figure('Position',[100 0 1000 1000])
    for c = 1:c_max
        disp(c)
        %Learning parameters
        del_W = zeros(N,N);
        del_M = zeros(N,N);
        ham = zeros(N,t_total,num_patterns);
        for p = randperm(num_patterns)
            [u_it, v_kt,rate_uit,rate_vkt,T, A_j] = deal(zeros(N,t_total));
            [del_ui] = deal(zeros(N,1));
            for t = 2:t_total/dt - (1/dt)
                sum_uj = W_ij*rate_uit(:,t-1);
                sum_ui = P_ik*rate_uit(:,t-1);
                sum_vk = M_ki*rate_vkt(:,t-1);
                max_inh = mean(rate_uit(:,t-1))-amp*f*N;
%                 del_ui = (chi(:,p,t) - u_it(:,t-1) + sum_uj)*(dt/tau_u);
                del_ui = (chi(:,p,t) - u_it(:,t-1) + sum_uj - sum_vk)*(dt/tau_u);
                del_vk = (chi(:,p,t) - v_kt(:,t-1) + sum_ui)*(dt/tau_v);
%                 del_vk = 0;
                u_it(:,t) = u_it(:,t-1) + del_ui;
                v_kt(:,t) = v_kt(:,t-1) + del_vk;
                for i = 1:N
                    rate_uit(i,t) = phi_2(u_it(i,t),v,theta_phi,u_c);
                    rate_vkt(i,t) = phi_2(v_kt(i,t),v,theta_phi,u_c);
                end
                H_t = eta_t*lambda2(rate_uit(:,t-1),thresh);
                del_T = (-T(:,t-1) + (H_t.*(T_max-T(:,t-1))))*(dt/tau_t);
                T(:,t) = T(:,t-1) + del_T;

                H_a = eta_a*chi(:,p,t-1);
    %             H_a = eta_a*rate_uit(:,t-1);
                del_A = (-A_j(:,t-1) + (H_a.*(A_max-A_j(:,t-1))))*(dt/tau_a);
                A_j(:,t) = A_j(:,t-1) + del_A;
    %             del_W = del_W + (eta_rec*T(:,t)*(A_j(:,t) - rate_uit(:,t).*A_t2(:,p,t))')*exp((-(t-t_reward(p)).^2)/20);
    %             del_W = del_W + (eta_rec*T(:,t)*(A_t2(:,p,t) - rate_uit(:,t).*A_t2(:,p,t))')*exp((-(t-t_reward(p)).^2)/20);
    %             del_W = del_W + (eta_rec*T(:,t)*(A_t2(:,p,t-1) - rate_uit(:,t-1).*(1-A_t2(:,p,t-1)))')*exp((-(t-100).^2)/20);
                del_W = del_W + (eta_rec*T(:,t)*lambda3((A_j(:,t) - rate_uit(:,t)),thresh3)')*exp((-(t-t_reward(p)).^2)/tau_exp);
                del_M = del_M + eta_inh*(lambda2(rate_vkt(:,t),thresh)*lambda2(rate_uit(:,t),thresh)');
                ham(:,t,p) = (A_t2(:,p,t) - rate_uit(:,t))*exp((-(t-t_reward(p)).^2)/tau_exp);
            end
            
            temp = abs(sum(ham(:,:,p),2));
            ham_ful(p,c) = 100*sum(temp)/(N*f*.5*sqrt(tau_exp*pi));
            W_ij = W_ij + del_W.*(abs(del_W)>thresh2);
            W_ij = W_ij.*(W_ij>0);
            W_ij = W_ij - diag(diag(W_ij));
            M_ki = M_ki - del_M + lambda*M_ki;
%             M_ki = 0;
            M_ki = M_ki.*(M_ki>0);
            M_ki(M_ki>1.5*f) = 1.5*f;
            rate_uit_ex{p} = rate_uit(A_t(:,p,t-1)>0,:);
            rate_uit_af(:,:,p) = rate_uit;
        end
        avg_rate = cellfun(@mean,rate_uit_ex,'UniformOutput',false);
        for i = 1:num_patterns
            subplot(3,num_patterns,i)
            imagesc(rate_uit_af(:,:,i))
            title(['HD is ', num2str(ham_ful(i,c))])
        end
% %     %     subplot(3,3,3)
% %     %     imagesc(ham(:,:,1))
% %     %     colorbar
% %     %     axis manual

        subplot(3,1,2)
        for i = 1:num_patterns
            plot(1:t_total,avg_rate{i})
            hold on
%             plot(1:t_total,T,'b')
%             hold on 
%             plot(1:t_total,A_j,'r*')
        end
% %     %     plot(1:t_total,10*mean(T,1))
% %     %     hold on 
% %     %     plot(1:t_total,10*mean(A_j,1))
        hold off
        ylim([0 15])
        subplot(3,2,5)
        imagesc(W_ij)
        title('W_ij')
        subplot(3,2,6)
        imagesc(M_ki)
        title('M_ki')
        drawnow

    end
    ph = 2*patterns-1;
    W_ij_hop = (1/N)*patterns*(patterns');
    crosstalk_hop = W_ij_hop*patterns - patterns;
    c_i_hop = -patterns.*crosstalk_hop;
    figure('Position',[100 0 1000 1000])
    %%
    k_length = 10;
    alpha = 10;
    hamming_pat = zeros(num_patterns,k_length);
    ham_ful2 = zeros(num_patterns,k_length);

    for k = 1:k_length

        chi_new = chi;
        for p = 1:num_patterns
    %         a = randperm(N,round((.2*f)*(N)));
            a = randperm(N,round(0.1*N));
            chi_new(a,p,1:100) = chi2(a,p,1:100);
        end
        rate_uit_ex2 = cell(1,num_patterns);
        rate_uit_af2 = zeros(N,t_total,num_patterns);
        [u_it, v_kt,rate_uit,rate_vkt,T, A_j] = deal(zeros(N,t_total));
        [del_ui] = deal(zeros(N,1));
        ham2 = zeros(N,t_total,num_patterns);
        for p = randperm(num_patterns)
            for t = 2:t_total/dt - (1/dt)
                sum_uj = W_ij*rate_uit(:,t-1);
                sum_ui = P_ik*rate_uit(:,t-1);
                sum_vk = M_ki*rate_vkt(:,t-1);
                if (t < 110) && (mean(rate_uit(:,t-1)) < amp*f)
                    max_inh = (amp*f + alpha) - mean(rate_uit(:,t-1));
                    max_inh = max_inh + (rand(50,1) - .5)/.1;
                elseif (t<110) && (mean(rate_uit(:,t-1)) > amp*f)
                    max_inh = (amp*f - alpha) - mean(rate_uit(:,t-1));
                    max_inh = max_inh + (rand(50,1) - .5)/.1;
                else
                    max_inh = 0;
                end
           
%                 max_inh = 0;
                del_ui = (chi_new(:,p,t) - u_it(:,t-1) + sum_uj + max_inh -sum_vk)*(dt/tau_u);
%                 del_ui = (chi_new(:,p,t) - u_it(:,t-1) + sum_uj + max_inh)*(dt/tau_u);
                del_vk = (chi_new(:,p,t) - v_kt(:,t-1)+sum_ui)*(dt/tau_v);
%                 del_vk = 0;
                u_it(:,t) = u_it(:,t-1) + del_ui;
                v_kt(:,t) = v_kt(:,t-1) + del_vk;
                for i = 1:N
                    rate_uit(i,t) = phi_2(u_it(i,t),v,theta_phi,u_c);
                    rate_vkt(i,t) = phi_2(v_kt(i,t),v,theta_phi,u_c);
                end
                ham2(:,t,p) = (A_t2(:,p,t) - rate_uit(:,t))*exp((-(t-t_reward(p)).^2)/tau_exp);
            end
            temp = abs(sum(ham2(:,:,p),2));
            ham_ful2(p,k) = 100*sum(temp)/(N*f*.5*sqrt(tau_exp*pi));
            rate_uit_ex2{p} = rate_uit(A_t(:,p,t-1)>0,:);
            rate_uit_af2(:,:,p) = rate_uit;
            hamming_pat_var = abs(rate_uit_af(:,:,p)-rate_uit_af2(:,:,p));
            hamming_pat(p,k) = mean(hamming_pat_var,'all');
        end
        avg_rate2 = cellfun(@mean,rate_uit_ex2,'UniformOutput',false);
        for i = 1:num_patterns
            subplot(2,num_patterns,i)
            imagesc(rate_uit_af2(:,:,i))
            title(['HD is ', num2str(ham_ful2(i,k))])
            HD_q((q-1)*num_patterns + i) = ham_ful2(i,k);
        end

        subplot(2,1,2)
        for i = 1:num_patterns
            plot(1:t_total,avg_rate2{i})
            hold on
        end
        hold off
        ylim([0 15])
        drawnow



    end

    total_ham = mean(hamming_pat,'all');
end

total_q = mean(HD_q<60);
disp(['Percentage of patterns recalled successfully is ', num2str(total_q*100), '%'])

function y = phi_2(u,v,theta,u_c)
if isnan(u)
    y = 0;
elseif theta >= u
    y = 0;
elseif (theta < u) && (u < u_c)
    y = v*(((u-theta)/(u_c-theta))^2);
elseif u_c <= u 
    y = 2*v*sqrt(((u-theta)/(u_c - theta))-.75);
end
end

function y = lambda2(u,thresh)
y = u.*(u>thresh);
end

function y = lambda3(u,thresh)
y = u.*(abs(u)>thresh);
end

function y = lambda4(u,thresh)
y = (u+thresh).*(u<=-thresh) + (u-thresh).*(u>=thresh);
end
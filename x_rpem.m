function varargout =...
    x_rpem(u,y,...
    nA,...
    nB,...
    strc_adapInit,...
    theta_init,...
    SW_TVgain)
% function varargout =...
%     max_nlms(u,w_d,...
%     adap_N,...
%     beta_init,beta_rate,beta_end,...
%     theta_init,...
%     adap_epsilon,...
%     SW_TVgain)

% u: noisy
% y: desired
% author: Xu Chen; mail.xchen@gmail.com
n_data = length(u);
if n_data~=length(y)
    error('Dimensions do not match')
end
%%
persistent N phi e_priori F y_priori y_post theta lambda thetaold phi_f yk_f uk_f;
% N       = uint16(nA+nB+1); % +1 because b_0 ~= 1
N       = uint16(nA+nB); % +1 because b_0 ~= 1
phi     = zeros(N,1);
phi_f     = zeros(N,1);
if size(theta_init,1)<size(theta_init,2)
    theta_init = theta_init';
end
theta   = theta_init;
F       = strc_adapInit.F_init;
y_priori   = 0;

% forgetting factor
lambda          = strc_adapInit.lambda_init;
if SW_TVgain == 1
    lambda_rate = strc_adapInit.lambda_rate;% 0.995;
    lambda_end  = strc_adapInit.lambda_end;% 0.99;
end
lambda_vec      = zeros(n_data,1);
%%
adap_para       = zeros(n_data,N);
tr_adap_gain    = zeros(1,n_data);
est_err = zeros(1,n_data);
ensemble_avg = zeros(1,n_data);
try
    window_L = strc_adapInit.window_L;
catch
    window_L = 100;
    disp('window_L not specified')
    p = mfilename('fullpath');
    edit([p,'.m'])
    keyboard
end
indx_sliding_window = 1;
cost_slidingWindow = zeros(1,n_data);
regres_sliding_window = zeros(window_L,1);
%%
for jj = 1:n_data
    %%
    % rpem has difficulty getting u(k) in phi_f(k) if the model is 
    %   b0 + b1 z^-1 + ...
    % the following was by assuming the model to be
    %   b1 z^-1 + ...
    % be careful about this when plotting the transfer function
    yk         = y(jj);
    uk         = u(jj);
    
    y_priori = phi'*theta;
    e_priori = yk - y_priori;
    %% performance evaluation
    regres_sliding_window(2:end) = regres_sliding_window(1:end-1);
    regres_sliding_window(1) = e_priori^2;
    
    if jj > 1
        ensemble_avg(jj) = 1/jj*( (jj-1)*ensemble_avg(jj-1) + e_priori^2 );
    else
        ensemble_avg(jj) = e_priori^2;
    end
    if 1
        if indx_sliding_window == window_L % 2012-04-11
            cost_slidingWindow(jj) = mean(regres_sliding_window);
            indx_sliding_window = 1;
        else
            if jj > 1
                cost_slidingWindow(jj) = cost_slidingWindow(jj-1);
            else
                cost_slidingWindow(jj) = 0;
            end
            %             if isempty(N_NORM)
            %                 N_NORM  = 40;
            %                 ITER    = 0;
            %                 u_vec   = zeros(N_NORM,1);
            %                 TMR_old = 0;
            %             end
            %             ITER            = ITER + 1;
            %             u_vec(2:N_NORM) = u_vec(1:N_NORM-1);
            %             u_vec(1)        = u;
            %             if ITER == N_NORM
            %                 TMR_old     = 3*sqrt(var(u_vec));
            %                 %     TMR_old     = 3*sqrt(mean(u_vec.^2));
            %                 ITER        = 0;
            %             end
            %             TMR             = TMR_old;
        end
        indx_sliding_window = indx_sliding_window + 1;
    else
        cost_slidingWindow(jj) = mean(regres_sliding_window);
    end
    % EOF performance evaluation
    %%
    lambda_vec(jj) = lambda;
    try
        if SW_TVgain == 1
            lambda = lambda_end - (lambda_end-lambda)*lambda_rate;
        elseif SW_TVgain == 2 % time dependedn only
        elseif SW_TVgain == 3 % slock's algorithm
        elseif SW_TVgain == 4
        elseif SW_TVgain == 5
        else
        end
    catch
        warning('Time-varying forgetting factor failed.')
    end
    F = 1/lambda*(F-F*phi_f*phi_f'*F/(lambda+phi_f'*F*phi_f));
    thetaold = theta;
    theta = theta + F*phi_f*e_priori/(lambda+phi_f'*F*phi_f); % PAA
    if any(abs(roots([1;theta(1:nA)]))>=1) % can be replaced by altorithms that are computationally more friendly
        theta = thetaold;
    end
    %% a posteriori
    y_post = theta'*phi;
    
    % phi: [-y(k-1),-y(k-2),...,-y(k-nA),...
    %       u(k),u(k-1),u(k-2),...,u(k-nB)]^T
    phi(2:nA)   = phi(1:nA-1);
    phi(1)      = -y_post;
%     phi(nA+2:nA+nB+1) = phi(nA+1:nA+nB);
    phi(nA+2:nA+nB) = phi(nA+1:nA+nB-1);
    phi(nA+1)   = uk;
    
%     yk_f = phi_yf'*theta(1:nA) + yk;
%     uk_f = phi_uf'*theta(1:nA) + uk;
    yk_f = phi_f(1:nA)'*theta(1:nA) + y_post;
    if nA == nB
        uk_f = -phi_f(nA+1:nA+nB)'*theta(1:nA) + uk;
    else
        error('tbd')
    end
    phi_f(2:nA) = phi_f(1:nA-1);
    phi_f(1) = -yk_f;
    phi_f(nA+2:nA+nB) = phi_f(nA+1:nA+nB-1);
    phi_f(nA+1)   = uk_f;
        
    est_err(jj) = e_priori;
    adap_para(jj,:) = theta(:)';
    tr_adap_gain(jj) = trace(F);
    if 0
        %%
        figure, plot(tr_adap_gain)
    end
end
%%
if nargout == 1
    varargout{1} = adap_para;
elseif nargout == 2
    varargout{1} = adap_para;
    varargout{2} = tr_adap_gain;
elseif nargout == 3
    varargout{1} = adap_para;
    varargout{2} = tr_adap_gain;
    varargout{3} = est_err;
    
elseif nargout == 4
    varargout{1} = adap_para;
    varargout{2} = tr_adap_gain;
    varargout{3} = est_err;
    varargout{4} = ensemble_avg;
    
elseif nargout == 5
    varargout{1} = adap_para;
    varargout{2} = tr_adap_gain;
    varargout{3} = est_err;
    varargout{4} = ensemble_avg;
    varargout{5} = cost_slidingWindow;
elseif nargout == 6
    varargout{1} = adap_para;
    varargout{2} = tr_adap_gain;
    varargout{3} = est_err;
    varargout{4} = ensemble_avg;
    varargout{5} = cost_slidingWindow;
    varargout{6} = lambda_vec;
elseif nargout == 0
    figure, plot(adap_para);
    ylabel 'parameter estimate'
    xlabel 'iteration'
    figure, plot(est_err)
    ylabel 'estimation error'
    xlabel 'iteration'
    figure, plot(tr_adap_gain)
    ylabel 'adaptation gain'
    xlabel 'iteration'
    figure, plot(ensemble_avg)
    ylabel 'ensemble average of the squared errors'
    xlabel 'iteration'
    figure, plot(lambda_vec)
    xlabel 'iteration'
    ylabel 'forgetting factor'
else
    error 'Error in the number of outputs.'
end
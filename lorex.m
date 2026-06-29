function [y_rec, IntermediateResults, FinalResults] = lorex(x, y, varargin)
% LOREX.m (Lorentzian Extractor) fits, and reconstructs multi-peak Lorentzian spectra.
%
% This function automates the detection and fitting of overlapped, 
% Lorentzian resonance peaks in 1D or 2D datasets.
%
% REQUIRED INPUTS:
%   x - 1D array of frequencies (Hz). Will be forced to column vector.
%   y - 1D array or 2D matrix of amplitudes. The number of rows in `y` 
%       must exactly match the length of `x`.
%
% OPTIONAL NAME-VALUE PAIRS:
%   'FWHM_guess'      - Estimated typical peak FWHM (Hz). (Default: 5*dx)
%   'AmpBounds'       - [min, max] allowed amplitude. (Default: [1e-6, max(y)*2])
%   'MinProminence'   - Fractional threshold for finding initial peaks. (Default: 0.01)
%   'DevTol'          - Allowed peak drift in Hz. (Default: 0.5)
%   'ResMult'         - Resolution multiplier for output arrays. (Default: 100)
%   'CheckPlots'      - Boolean to show step-by-step diagnostic plots. (Default: false)
%   'ShowHeatmap'     - Boolean to show final 2D before/after heatmaps. (Default: true)
%
% OUTPUTS:
%   y_rec               - The reconstructed, high-res. fitted data.
%   IntermediateResults - MATLAB table of the initial guesses.
%   FinalResults        - MATLAB table of the final multi-peak fit.
%
% -------------------------------------------------------------------------
% Author: [Gengming Liu]
% Date: June 2026
% Institution: [University of Illinois Urbana-Champaign]
% Contact: [hi@gengming.me]
% Version: 1.0.0
%
% Citation: 
% If you use this code in your research, please cite it via Zenodo:
% https://doi.org/10.5281/zenodo.21019679
%
% Acknowledgments:
% This code includes the 'lorentzfit.m' algorithm by Jered R. Wells (BSD License)
% retrieved from the MATLAB Central File Exchange.
% https://www.mathworks.com/matlabcentral/fileexchange/33775-lorentzfit-x-y-varargin
% -------------------------------------------------------------------------

    %% 1. Parse Inputs and Enforce Dimensions
    p = inputParser;
    
    addRequired(p, 'x', @(a) isnumeric(a) && isvector(a));
    addRequired(p, 'y', @isnumeric);
    
    addParameter(p, 'FWHM_guess', [], @isnumeric);
    addParameter(p, 'AmpBounds', [], @(a) isnumeric(a) && length(a)==2);
    addParameter(p, 'MinProminence', 0.01, @isnumeric);
    addParameter(p, 'DevTol', 0.5, @isnumeric);
    addParameter(p, 'ResMult', 100, @isnumeric);
    addParameter(p, 'CheckPlots', false, @islogical);
    addParameter(p, 'ShowHeatmap', true, @islogical);
    
    parse(p, x, y, varargin{:});
    opts = p.Results;
    
    % Enforce X as column vector
    x = x(:);
    Nx = length(x);
    
    % Check for uniform spacing (Required for 5-point stencil derivative)
    dx_array = diff(x);
    if (max(dx_array) - min(dx_array)) > 1e-4 * mean(dx_array)
        error('The frequency vector (x) must be uniformly spaced. Please interpolate your data onto a uniform grid before fitting.');
    end
    dx = mean(dx_array);
    
    % Enforce Y orientation (Strict)
    if isvector(y)
        y = y(:);
    end
    
    if size(y, 1) ~= Nx
        error('Dimension mismatch: The number of rows in y (%d) must match the length of x (%d). If your data is transposed, please transpose y before passing it to this function.', size(y,1), Nx);
    end
    num_data_traces = size(y, 2);

    %% 2. Establish Dynamic Defaults & Bounds
    if isempty(opts.FWHM_guess)
        opts.FWHM_guess = 5 * dx; 
    end
    if isempty(opts.AmpBounds)
        opts.AmpBounds = [1e-6, max(y(:)) * 2]; 
    end
    
    % Auto-calculate bounds based on the FWHM guess
    FWHM_bounds = [dx, 2 * opts.FWHM_guess];
    gamma_bounds = (FWHM_bounds ./ 2).^2;
    
    % Strictly determine DataPoints window based on guess and step size
    DataPoints = round(opts.FWHM_guess / dx);
    if mod(DataPoints, 2) == 0
        DataPoints = DataPoints + 1; 
    end
    DataPoints = max(3, DataPoints - 2); 
    halfw = (DataPoints - 1) / 2;
    
    % Strict silent options for the individual lorentzfit calls
    silent_opts = optimset('Display', 'off');

    %% 3. Setup Outputs
    finef = linspace(x(1), x(end), Nx * opts.ResMult);
    y_rec = zeros(size(y));
    
    rows1 = struct('data_trace',{},'peak',{},'x0',{},'gamma',{},'amplitude',{});
    rows2 = struct('data_trace',{},'peak',{},'x0',{},'gamma',{},'amplitude',{});

    %% 4. Core Fitting Loop
    for col = 1:num_data_traces
        chain = y(:, col);
        y2 = nan(length(chain),1);  

        % 5-point stencil second derivative
        n = 3:(length(chain)-2);
        y2(n) = (-chain(n+2) + 16*chain(n+1) - 30*chain(n) + 16*chain(n-1) - chain(n-2)) ./ (12*dx^2);
        d2 = -1 .* y2(3:(end-2)); 

        [pks1, locs1, w1, ~] = findpeaks(d2, 'MinPeakProminence', opts.MinProminence * max(d2));

        area = zeros(length(locs1),1);
        d22 = max(d2, 0);

        % Area calculation with dynamic window and boundary safety
        for xx = 1:length(locs1)
            idx_center = locs1(xx);
            idx_start = max(1, idx_center - halfw);
            idx_end   = min(length(d22), idx_center + halfw);
            
            x_win = x(idx_start+2 : idx_end+2);
            y_win = d22(idx_start : idx_end);
            area(xx) = trapz(x_win, y_win);
        end

        dumppeak = area > (max(area) * opts.MinProminence);
        
        pks = nonzeros(pks1 .* dumppeak);
        locs = nonzeros(locs1 .* dumppeak);
        w = nonzeros(w1 .* dumppeak);

        if opts.CheckPlots
            fig_data = figure('Name', sprintf('Diagnostics: Data Trace %d', col));
            
            % --- Subplot 1: Raw Data & Peak Finding ---
            subplot(3,1,1)
            h_raw = plot(x, chain(:), 'k.'); hold on;
            xlabel('Frequency'); ylabel('Amplitude');
            title(sprintf('Measured #%d', col));
            ylim([-0.1 max(chain)*1.1]); xlim([x(1) x(end)]);
            
            for xv = x(locs1+2)', xline(xv, 'g--'); end
            for xv = x(locs+2)',  xline(xv, 'r--'); end
            
            % Dummy plots for robust legends
            h_cand = plot(NaN, NaN, 'g--');
            h_kept = plot(NaN, NaN, 'r--');
            legend([h_raw, h_cand, h_kept], {'Raw Data', 'Dropped Peaks', 'Kept Peaks'});

            % --- Separate Figure: d2y Area Analysis (No red lines) ---
            figure('Name', sprintf('Peak Analysis: Data Trace %d', col));
            h_d2 = plot(x(3:end-2), d22, 'b-'); hold on;
            h_area = plot(x(locs1+2), area, 'ok', 'MarkerFaceColor', 'k');
            h_thresh = yline(max(area) * opts.MinProminence, 'c--', 'LineWidth', 1.5);
            
            xlabel('Frequency'); ylabel('Second Derivative');
            title(sprintf('Trace #%d Peak Rejection Analysis', col));
            legend([h_d2, h_area, h_thresh], {'d2y (2nd Deriv)', 'Candidate Areas', 'Area Threshold'});
        end

        % ---- Initial individual fits ----
        leftlocs = zeros(length(pks),1);
        for i = 1:length(pks)
            iiL = max(1, locs(i) + 2 - halfw);
            iiR = min(Nx, locs(i) + 2 + halfw);
            tofitx = x(iiL:iiR);
            tofity = chain(iiL:iiR);

            gamma_g = w(i)/10;
            p0 = [pks(i), x(locs(i)+2), (gamma_g/2).^2];
            
            lb = [opts.AmpBounds(1), x(1), gamma_bounds(1)];
            ub = [opts.AmpBounds(2), x(end), gamma_bounds(2)];
            
            [~, params] = lorentzfit(tofitx, tofity, p0, [lb; ub], '3', silent_opts);

            isWeak = (params(1) < opts.AmpBounds(1));

            if ~isWeak 
                leftlocs(i) = 1;
                rows1(end+1) = struct('data_trace', col, 'peak', i, ...
                                      'x0', params(2), 'gamma', params(3), 'amplitude', params(1));
            end
        end  
        w = nonzeros(w .* leftlocs);

        % ---- Refit deviated individual peaks on residuals ----
        thisResRows = [rows1.data_trace] == col;
        idx_rows = find(thisResRows);
        r1 = rows1(idx_rows);

        A0  = [r1.amplitude]'; x00 = [r1.x0]'; g0  = [r1.gamma]';
        K = numel(A0);

        if K > 0
            x00_guess = nonzeros(x(locs + 2) .* leftlocs);
            dev_mask = abs(x00 - x00_guess) > opts.DevTol;
            dev_idx  = find(dev_mask).';
            good_idx = find(~dev_mask).';

            xx = x(:);
            y_good = zeros(numel(xx),1);
            if ~isempty(good_idx)
                den_g  = (xx - x00(good_idx).').^2 + g0(good_idx).'; 
                y_good = sum(bsxfun(@rdivide, A0(good_idx).', den_g), 2);
            end

            y_refit = zeros(numel(xx),1);

            for k = dev_idx
                resid_full = chain(:) - y_good - y_refit;
                iiL = max(locs(k)+2-halfw, 1);
                iiR = min(locs(k)+2+halfw, numel(x));
                tofitx2 = x(iiL:iiR);
                tofity2 = max(resid_full(iiL:iiR), 0);
                tofity2(~isfinite(tofity2)) = 0;

                gamma_g = w(k)/10;                                   
                p0k = [max(tofity2), x00_guess(k), (gamma_g/2).^2];
                
                lb = [opts.AmpBounds(1), x(1), gamma_bounds(1)];
                ub = [opts.AmpBounds(2), x(end), gamma_bounds(2)];
                
                try
                    [~, par2] = lorentzfit(tofitx2, tofity2, p0k, [lb; ub], '3', silent_opts);
                    A0(k)  = max(par2(1), eps);
                    x00(k) = par2(2);
                    g0(k)  = max(par2(3), eps);

                    rows1(idx_rows(k)).amplitude = A0(k);
                    rows1(idx_rows(k)).x0        = x00(k);
                    rows1(idx_rows(k)).gamma     = g0(k);
                    
                    y_refit = y_refit + A0(k) ./ ((xx - x00(k)).^2 + g0(k));
                catch
                end
            end
        end

        % ---- Summed constrained fit (Final Pass) ----
        [A, x0, g, ~] = fit_multi_lorentz_dump(x, chain, A0, x00, g0, opts.DevTol, opts.AmpBounds, gamma_bounds);

        for i = 1:length(A)
           rows2(end+1) = struct('data_trace', col, 'peak', i, 'x0', x0(i), 'gamma', g(i), 'amplitude', A(i)); 
        end
        
        if opts.CheckPlots
            % --- Subplot 2: Individual Fits ---
            figure(fig_data); subplot(3,1,2); hold on;
            h_raw2 = plot(x, chain(:), 'k.');
            
            den_final = (finef(:) - x0(:).').^2 + g(:).';
            Ypk_final = bsxfun(@rdivide, A(:).', den_final); 
            
            Kfin = size(Ypk_final,2);
            skip = max(1, floor(numel(finef)/4000));
            idx  = 1:skip:numel(finef);
            Xc   = [repmat(finef(idx).', 1, Kfin); nan(1,Kfin)];
            Yc   = [Ypk_final(idx,:); nan(1,Kfin)];
            h_indiv = plot(Xc(:), Yc(:), 'c--');                     
            
            for xv = x0', xline(xv, 'r-'); end
            h_center = plot(NaN, NaN, 'r-'); % Dummy for legend
            
            xlim([x(1) x(end)]); ylim([-0.1 max(chain)*1.1]);
            xlabel('Frequency'); ylabel('Amplitude'); title('Individual Peaks');
            legend([h_raw2, h_indiv, h_center], {'Raw Data', 'Individual Peaks', 'Peak Centers'});

            % --- Subplot 3: Final Summation ---
            subplot(3,1,3); hold on;
            h_raw3 = plot(x, chain(:), 'bo');
            h_sum = plot(finef, sum(Ypk_final,2), 'r-');
            
            xlim([x(1) x(end)]); ylim([-0.1 max(chain)*1.1]);
            xlabel('Frequency'); ylabel('Amplitude'); title('Overall Comparision');
            legend([h_raw3, h_sum], {'Raw Data', 'Fits Summed'});
        end
        
        y_rec(:, col) = sum( bsxfun(@rdivide, A(:).', (x(:) - x0(:).').^2 + g(:).'), 2 );
    end

    IntermediateResults = struct2table(rows1);
    FinalResults = struct2table(rows2);

    %% 5. Plot Side-by-Side Heatmaps
    if opts.ShowHeatmap
        figure('Name', 'Heatmap Comparison', 'Position', [100, 100, 1000, 400]);
        
        % Original Heatmap
        subplot(1,2,1);
        imagesc((y)); shading flat; colorbar; colormap(jet);
        set(gca,'Ydir','normal');
        set(gca,'YTick',1:10:Nx);
        set(gca,'YTickLabel', round(x(1:10:end),1));
        xlabel('Data Dimension'); ylabel('Frequency');
        title('Original');
        
        % Reconstructed Heatmap
        subplot(1,2,2);
        imagesc((y_rec)); shading flat; colorbar; colormap(jet);
        set(gca,'Ydir','normal');
        set(gca,'YTick',1:10:Nx);
        set(gca,'YTickLabel', round(x(1:10:end),1));
        xlabel('Data Dimension'); ylabel('Frequency');
        title('Fitted');
    end
end

%% =====================================================================
%  LOCAL HELPER FUNCTIONS
%  =====================================================================

function [A,x0,g,yfit] = fit_multi_lorentz_dump(x, y, A0, x00, g0, DevTol, AmpBounds, gamma_bounds)
    % FIT_MULTI_LORENTZ_DUMP fits a sum of Lorentzians, dropping degenerate peaks dynamically.
    MaxRefits = 2;
    AmpTolRel = 1e-4;   
    
    opts_lsq = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxFunctionEvaluations', 5e4, 'MaxIterations', 2e3, 'FunctionTolerance', 1e-12);

    x = x(:); y = y(:);
    A0 = A0(:); x00 = x00(:); g0 = g0(:);
    Kall = numel(A0);
    m = numel(y);
    
    lorentz_sum = @(xx,A,x0,g) sum(A.' ./ ((xx - x0.').^2 + g.'), 2);
    model  = @(pp,xx,K) lorentz_sum(xx, pp(1:K), pp(K+(1:K)), pp(2*K+(1:K)));
    pack   = @(A,x0,g) [A; x0; g];
    unpack = @(pp,K) deal(pp(1:K), pp(K+(1:K)), pp(2*K+(1:K)));

    kept_mask = true(Kall,1);
    dropped_mask = false(Kall,1);

    for refit_iter = 1:MaxRefits
        keep_idx = find(kept_mask);
        K = numel(keep_idx);
        if K==0
            A=[]; x0=[]; g=[]; yfit=zeros(size(x)); return;
        end
        if m < 3*K
            error('Not enough data points: m=%d, need >=3K=%d.', m, 3*K);
        end

        A0k  = A0(keep_idx); x0k0 = x00(keep_idx); g0k = g0(keep_idx);
        dx   = median(diff(x));

        lbA = AmpBounds(1) * ones(K,1);
        ubA = AmpBounds(2) * ones(K,1);
        lbg = gamma_bounds(1) * ones(K,1);
        ubg = gamma_bounds(2) * ones(K,1);

        slack = max(DevTol, 2*dx);     
        lbx = x0k0 - slack;
        ubx = x0k0 + slack;

        lb = [lbA; lbx; lbg];
        ub = [ubA; ubx; ubg];

        g0k  = min(max(g0k, lbg + 0.25*(ubg-lbg)), lbg + 0.75*(ubg-lbg));
        p0   = pack(A0k,x0k0,g0k);

        modelK = @(pp,xx) model(pp,xx,K);

        [p_hat_k,~,~,~] = lsqcurvefit(modelK,p0,x,y,lb,ub,opts_lsq);
        [Ak,x0k_hat,gk] = unpack(p_hat_k,K);

        dev = abs(x0k_hat - x0k0) > DevTol;
        amp_low = Ak <= (max(A0k)*AmpTolRel);
        offenders_local = find(dev | amp_low);

        if isempty(offenders_local)
            A=Ak; x0=x0k_hat; g=gk; yfit=modelK(p_hat_k,x); return;
        else
            offenders_global = keep_idx(offenders_local);
            kept_mask(offenders_global)    = false;
            dropped_mask(offenders_global) = true;
        end
    end

    keep_idx = find(kept_mask);
    K = numel(keep_idx);
    if K==0
        A=[]; x0=[]; g=[]; yfit=zeros(size(x)); return;
    end

    A0k  = A0(keep_idx); x0k0 = x00(keep_idx); g0k = g0(keep_idx);
    
    lbA = AmpBounds(1) * ones(K,1);
    ubA = AmpBounds(2) * ones(K,1);
    lbg = gamma_bounds(1) * ones(K,1);
    ubg = gamma_bounds(2) * ones(K,1);
    lbx = x0k0 - max(DevTol, 2*median(diff(x)));
    ubx = x0k0 + max(DevTol, 2*median(diff(x)));

    lb = [lbA; lbx; lbg];
    ub = [ubA; ubx; ubg];

    g0k = min(max(g0k, lbg + 0.25*(ubg-lbg)), lbg + 0.75*(ubg-lbg));
    p0  = pack(A0k,x0k0,g0k);

    modelK = @(pp,xx) model(pp,xx,K);
    [p_hat,~,~,~] = lsqcurvefit(modelK,p0,x,y,lb,ub,opts_lsq);
    [A,x0,g] = unpack(p_hat,K);

    amp_low = A <= (max(A0k)*AmpTolRel);
    if any(amp_low)
        A = A(~amp_low); x0 = x0(~amp_low); g = g(~amp_low);
    end
    yfit = lorentz_sum(x, A, x0, g);
end

%% =====================================================================
% LORENTZFIT BY JERED R. WELLS (BSD LICENSE INTACT)
% =====================================================================
function varargout = lorentzfit(x,y,varargin)
% Jered R Wells
% 11/15/11
% jered [dot] wells [at] duke [dot] edu
% v1.7 (2020/06/15)
narginchk(2,6); nargoutchk(0,5); fname = 'lorentzfit';
inputcheck(x,{'numeric'},{'real','nonnan','nonempty','finite'},fname,'X',1);
inputcheck(y,{'numeric'},{'real','nonnan','nonempty','finite','size',size(x)},fname,'Y',2);
p3 = ((max(x(:))-min(x(:)))./10).^2; p2 = (max(x(:))+min(x(:)))./2; p1 = max(y(:)).*p3; c = min(y(:));
optargs = {[],[],'3c',optimset('TolFun',max(mean(y(:))*1e-6,1e-15),'TolX',max(mean(x(:))*1e-6,1e-15),'Display','off')};
numvarargs = length(varargin);
for ii = 1:numvarargs; if isempty(varargin{ii}); varargin{ii} = optargs{ii}; end; end
optargs(1:numvarargs) = varargin;
[p0,bounds,nparams,options] = optargs{:};
if ~isempty(p0); inputcheck(p0,{'numeric'},{'real','nonnan','vector'},fname,'P0',3); end
if ~isempty(bounds)
    inputcheck(bounds,{'numeric'},{'real','nonnan','nrows',2},fname,'BOUNDS',4);
    lb = bounds(1,:); ub = bounds(2,:);
else
    lb = []; ub = [];
end
inputcheck(nparams,{'char'},{},fname,'NPARAMS',5); inputcheck(options,{'struct'},{},fname,'OPTIONS',6);
switch lower(nparams)
    case '3'
        if isempty(p0); p0 = [p1 p2 p3]; end
        [params,resnorm,residual,~,~,~,J] = lsqcurvefit(@lfun3,p0,x,y,lb,ub,options);
        yprime = lfun3(params,x);
    otherwise
        error('Custom implementation truncated for space: Only model 3 supported in this wrapper.');
end
varargout = {yprime,params,resnorm,residual,J};
end
function F = lfun3(p,x)
F = p(1)./((x-p(2)).^2+p(3));
end

function varargout = inputcheck(A,varargin)
V = true; ME = MException('','');
classes = varargin{1}; attributes = varargin{2};
try validateattributes(A,classes,attributes); catch ME; V = false; end
varargout = {V,ME};
end

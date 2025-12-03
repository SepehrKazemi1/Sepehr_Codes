%% CRU vs Synoptic (IRIMO) — RMSE & R^2 per station, subplots by region
clc; clear; close all;

fileCRU  = ""; %monthly cru output
fileIRIMO= ""; %use ur synoptic monthly averaged parameters for tmax and tmin and summed up prc and pet
faStationsOpt = "";   % u can rename ur stations here , if u dont want to, just make an excel with same names
outDir = "";  % 
    
varSheets = {'prc','pet','tmax','tmin'};

% choose the region u want ur stations within
Top = 33.696098; Bottom = 31.195833; Left = 50.033333; Right = 53.400000;
midLat = (Top+Bottom)/2; midLon = (Left+Right)/2;

% change font and colors
set(groot,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
clrSyn = [0.00 0.45 0.74]; % blue
clrCRU = [0.85 0.33 0.10]; % orange
regColor = struct('NW',[0.00 0.45 0.74], 'NE',[0.85 0.33 0.10], ...
                  'SW',[0.47 0.67 0.19], 'SE',[0.85 0.10 0.10]);

unitTxt = struct('prc','mm/mon','pet','mm/mon','tmax','°C','tmin','°C');

%for farsi stations
haveFA = false; S = table();
if isfile(faStationsOpt)
    try
        S = readtable(faStationsOpt,'PreserveVariableNames',true);
        haveFA = true;
    catch, haveFA = false;
    end
end
getFA  = @(en) get_fa_name(S, en);
getLat = @(en) get_coord(S, en, 'lat');
getLon = @(en) get_coord(S, en, 'lon');
LRM = char(8206);   % left-to-right Mark
RLM = char(8207);   %right-to-left Mark

for v = 1:numel(varSheets)
    varName = varSheets{v};

    Tc = readtable(fileCRU,   'Sheet', varName, 'VariableNamingRule','preserve');
    Ti = readtable(fileIRIMO, 'Sheet', varName, 'VariableNamingRule','preserve');

    tc = normalize_month_time(Tc{:,1});
    ti = normalize_month_time(Ti{:,1});

    namesC = string(Tc.Properties.VariableNames(2:end));
    namesI = string(Ti.Properties.VariableNames(2:end));

    % header mapping 
    map = struct('en',{},'ic',{},'ii',{},'fa',{});
    for k = 1:numel(namesI)
        en = namesI(k);
        j  = find(name_norm(namesC) == name_norm(en), 1);
        if ~isempty(j)
            rec.en = en; rec.ic = j+1; rec.ii = k+1;  % +1 for date column
            rec.fa = getFA(en);
            map(end+1) = rec; %#ok<SAGROW>
        end
    end
    if isempty(map), warning('No matching station headers for "%s".', varName); continue; end

    % for regionalization
    regions = {'NW','NE','SW','SE'};
    ByReg = struct('NW',[],'NE',[],'SW',[],'SE',[]);
    for iMap = 1:numel(map)
        en = map(iMap).en;
        if haveFA
            lat = getLat(en); lon = getLon(en);
            if ~isfinite(lat) || ~isfinite(lon), continue; end
            if ~(lon>=Left && lon<=Right && lat>=Bottom && lat<=Top), continue; end
            if     lat>=midLat && lon<=midLon, reg='NW';
            elseif lat>=midLat && lon> midLon, reg='NE';
            elseif lat< midLat && lon<=midLon, reg='SW';
            else,                                reg='SE';
            end
        else
            reg = 'NW';
        end
        ByReg.(reg) = [ByReg.(reg), map(iMap)]; %#ok<AGROW>
    end

    % plot per region 
    for r = 1:numel(regions)
        reg = regions{r};
        entries = ByReg.(reg);
        if isempty(entries), continue; end

        N = numel(entries);
        nrow = ceil(sqrt(N)); ncol = ceil(N/nrow);

        fig = figure('Color','w','Position',[80 80 1280 820], ...
            'Name', sprintf('%s — %s (CRU vs Syn)', reg, upper(varName)));
        tl  = tiledlayout(nrow,ncol,'Padding','compact','TileSpacing','compact');

        h1 = gobjects(0); h2 = gobjects(0); % for legend

        for k = 1:N
            ax = nexttile(tl, k);
            en = entries(k).en;
            fa = entries(k).fa;  nameToShow = ifelse(strlength(fa)>0, fa, en);

            % align by time intersection
            yC = double(Tc{:, entries(k).ic});
            yI = double(Ti{:, entries(k).ii});
            [t, iC, iI] = intersect(tc, ti);
            yC = yC(iC); yI = yI(iI);

            % keep only months present in Synoptic + finite on both
            m = isfinite(yI) & isfinite(yC);
            t = t(m); yC = yC(m); yI = yI(m);

            % metrics (RMSE + regression R^2)
            [rmse, r2] = metrics_reg(yI, yC);

            % robust axis limits
            allVals = [yI; yC];
            fin = isfinite(allVals);
            if any(strcmpi(varName,{'prc','pet'}))
                ylo = 0;
                if any(fin)
                    yhi = prctile(allVals(fin), 99);
                    if ~isfinite(yhi) || yhi<=0, yhi = max(allVals(fin)); end
                else
                    yhi = 1;
                end
                yhi = 1.2*yhi;
            else % tmax/tmin
                if any(fin)
                    q1 = prctile(allVals(fin),1);
                    q9 = prctile(allVals(fin),99);
                    span = q9 - q1;
                    if ~isfinite(span) || span<=0
                        mval = mean(allVals(fin)); sd = std(allVals(fin)); if sd==0, sd=1; end
                        ylo = mval - 3*sd; yhi = mval + 3*sd;
                    else
                        ylo = q1 - 0.1*span; yhi = q9 + 0.1*span;
                    end
                else
                    ylo = -1; yhi = 1;
                end
            end

            % plot
            hold(ax,'on');
            p1 = plot(ax, t, yI, '-o', 'Color', clrSyn, 'LineWidth',1.0, 'MarkerSize',3, 'DisplayName','Synoptic');
            p2 = plot(ax, t, yC, '-s', 'Color', clrCRU, 'LineWidth',1.0, 'MarkerSize',3, 'DisplayName','CRU');
            grid(ax,'on'); box(ax,'on');
            set(ax,'YLim',[ylo yhi], 'XLim',[min(t) max(t)]);
            ylabel(ax, unitTxt.(varName));
            metricsStr = sprintf('%sRMSE=%.2f   R^{2}=%.3f%s', LRM, rmse, r2, LRM);  % force LTR
            nameStr    = [RLM char(nameToShow)];                                      % force RTL
            title(ax, [metricsStr '  |  ' nameStr], ...
                  'Interpreter','tex', ...
                  'Color', regColor.(reg), 'FontWeight','bold');


            if k==1, h1=p1; h2=p2; end
        end

        xlabel(tl, 'Year (monthly)');
        legend([h1 h2], {'Synoptic','CRU'}, 'NumColumns',2, 'Location','southoutside');
        title(tl, sprintf('%s — %s  (CRU vs Synoptic).  n=%d', reg, upper(varName), N), ...
              'Color', regColor.(reg), 'FontWeight','bold');
        drawnow;  % make sure everything is rendered

        % Safe filename: e.g., "PRC_NW.png"
        fname = sprintf('%s_%s.png', upper(varName), reg);
        outPath = fullfile(outDir, fname);
        
        % Best with tiledlayout (R2020a+):
        try
            exportgraphics(fig, outPath, 'Resolution',300, 'BackgroundColor','white');
        catch
            % Fallback for older MATLAB
            set(fig,'PaperPositionMode','auto');
            print(fig, outPath, '-dpng', '-r300');
        end
        
        % Optional: close the figure to avoid a pile-up
        % close(fig)

    end
end

%% funcs
function [rmse, R2] = metrics_reg(y, x)
    % RMSE
    rmse = sqrt(mean((y - x).^2, 'omitnan'));
    % regression R^2 for y ~ a + b*x (with intercept)
    mask = isfinite(y) & isfinite(x);
    y = y(mask); x = x(mask);
    if numel(y) < 2, R2 = NaN; return; end
    X = [ones(size(x)) x];
    b = X\y;                % OLS fit
    yhat = X*b;
    sst = sum((y - mean(y)).^2);
    ssr = sum((yhat - mean(y)).^2);
    R2  = ssr / max(sst, eps);
end

function t = normalize_month_time(v)
    if isdatetime(v), t = dateshift(v,'start','month'); return; end
    if isnumeric(v)
        try, t = datetime(v,'ConvertFrom','excel'); catch, t = datetime(v,'ConvertFrom','posixtime'); end
        t = dateshift(t,'start','month'); return;
    end
    vs = string(v); ok=false;
    fmts = {'yyyy-MM','yyyy/MM','dd/MM/yyyy','MM/dd/yyyy','MMM d, yyyy','MMMM d, yyyy'};
    for f = fmts
        try, t = datetime(vs,'InputFormat',f{1}); ok=true; break; end %#ok<CCAT>
    end
    if ~ok, t = datetime(vs); end
    t = dateshift(t,'start','month');
end

function n = name_norm(x)
    x = lower(string(x));
    n = regexprep(x, '[\s_\-\(\)‌]', '');
end

function s = ifelse(c,a,b), if c, s=a; else, s=b; end, end

function fa = get_fa_name(S, en)
    fa = "";
    if isempty(S), return; end
    if all(ismember({'station_name','station_name_fa'}, string(S.Properties.VariableNames)))
        idx = find(strcmpi(string(S.station_name), string(en)), 1);
        if ~isempty(idx), fa = string(S.station_name_fa(idx)); end
    end
end
function val = get_coord(S, en, which)
    val = NaN;
    if isempty(S), return; end
    if any(strcmpi(string(S.Properties.VariableNames), 'station_name')) && ...
       any(strcmpi(string(S.Properties.VariableNames), which))
        idx = find(strcmpi(string(S.station_name), string(en)), 1);
        if ~isempty(idx), val = S.(which)(idx); end
    end
end

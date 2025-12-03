%%CRU monthly pdfs by quadrant (prc/pet/tmin/tmax) values in mm/month for prc and pet
clc; clear; close all;

dirCRU = "C:\Users\Asus\Desktop\CRU_REFINED";
outDir = "";
varsUser = {'prc','pet','tmin','tmax'};
map2CRU  = struct('prc','pre','pet','pet','tmin','tmn','tmax','tmx');

xLabelFA = struct( ...
  'prc' , 'بارش (میلی‌متر/ماه)','pet' , 'تبخیر-تعرق مرجع (میلی‌متر/ماه)','tmin', 'دمای کمینه (°C)','tmax', 'دمای بیشینه (°C)' );
yLabelFA = 'چگالی احتمال';

%region colors
regColor = struct('NW',[0.00 0.45 0.74], 'NE',[0.85 0.33 0.10], 'SW',[0.47 0.67 0.19],  'SE',[0.85 0.10 0.10]);

%bounds
Top = 33.696098; Bottom =31.195833; Left =50.033333; Right =53.400000;
midLat = (Top+Bottom)/2; midLon =(Left+Right)/2;

% period
Y0 =1980; Y1 =2024;    %monthly CRU period
thrPrc =1;              %mm/month threshold for "dry" vs "wet"

%figurenots 
set(groot,'DefaultAxesFontName','Arial','DefaultAxesFontWeight','bold');
set(groot,'DefaultTextFontName','Arial');

%main loop
for v =1:numel(varsUser)
    vs      =varsUser{v};
    vCRU    =map2CRU.(vs);
    ncPath  =findCruFile(dirCRU, vCRU);

    % read coords & time
    lon  =ncread(ncPath,'lon');
    lat  =ncread(ncPath,'lat');
    tnum =ncread(ncPath,'time');
    tun  =ncreadatt(ncPath,'time','units');
    time = cru_time_to_datetime(tnum, tun);   %start month

    %long
    if max(lon) > 180
        lon = mod(lon+180,360)-180; [lon, iOrd] =sort(lon);
    else
        iOrd = 1:numel(lon);
    end

    % spatial subset
    ilon = find(lon >= Left & lon <= Right);
    ilat = find(lat >= Bottom & lat <= Top);
    % time window subset
    maskT = (time >= datetime(Y0,1,1)) &(time <= datetime(Y1,12,1));
    tStart = find(maskT,1,'first');
    tEnd   = find(maskT,1,'last');
    nT     = tEnd - tStart + 1;

    info = ncinfo(ncPath, vCRU);
    dimNames = string({info.Dimensions.Name});
    [start,count] = deal(ones(1,numel(dimNames)));
    for d = 1:numel(dimNames)
        switch lower(dimNames(d))
            case 'time', start(d)=tStart;            count(d)=nT;
            case 'lat',  start(d)=min(ilat);         count(d)=numel(ilat);
            case 'lon',  start(d)=min(iOrd(ilon));   count(d)=numel(ilon);
        end
    end
    V = ncread(ncPath, vCRU, start, count);  %raw units
    ord = lower(dimNames);
    [~,iT]  = ismember("time",ord);
    [~,iLa] = ismember("lat",ord);
    [~,iLo] = ismember("lon",ord);
    V = permute(V, [iT iLa iLo]);

    lonSub = lon(ilon); latSub = lat(ilat);
    if max(ncread(ncPath,'lon')) > 180
        [lonSub, sortIdx] = sort(lonSub);
        V = V(:,:,sortIdx);
    end

    %unit change
    if strcmp(vs,'pet')
        tsub = time(tStart:tEnd);
        dpm  = eomday(year(tsub), month(tsub));     % day per monts
        scale = reshape(dpm, [numel(tsub) 1 1]);    % broadcast to [T x 1 x 1]
        V = V .* scale;                              % mmday to mm/month
    end

%4 region split
    [Lon2,Lat2] = meshgrid(lonSub, latSub);
    maskNW = (Lat2 >= midLat) & (Lon2 <= midLon);
    maskNE = (Lat2 >= midLat) & (Lon2 >  midLon);
    maskSW = (Lat2 <  midLat) & (Lon2 <= midLon);
    maskSE = (Lat2 <  midLat) & (Lon2 >  midLon);
    regs   = {'NW','NE','SW','SE'}; masks = {maskNW, maskNE, maskSW, maskSE};
%collect values 
    dataReg = struct('NW',[],'NE',[],'SW',[],'SE',[]);
    tmp = reshape(V, nT, []); % T x (Y*X)
    for r = 1:4
        mk = masks{r};
        if ~any(mk,'all'), dataReg.(regs{r}) = []; continue;
        end
        y = tmp(:, find(mk)); y = y(:);
        y = y(isfinite(y));
        if ismember(vs, {'prc','pet'}), y = y(y>=0);
        end
        dataReg.(regs{r}) = y;
    end

    % ---------------- Plot ----------------
    fig = figure('Color','w','Position',[80 80 1100 820], ...
                 'Name', sprintf('PDF %s  %d-%d', upper(vs), Y0, Y1));
    tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    for r = 1:4
        ax = nexttile(tl, r); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        y = dataReg.(regs{r});
        if numel(y) < 5
            text(ax,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center'); axis(ax,'off'); continue;
        end
        [xi,fi] = compute_pdf(y); plot(ax, xi, fi, 'LineWidth', 1.8);
        xlo = prctile(xi,1); xhi = prctile(xi,99);
        if ~isfinite(xlo)||~isfinite(xhi)||xlo==xhi, xlo=min(xi); xhi=max(xi); end
        xlim(ax,[xlo xhi]); ylabel(ax, yLabelFA, 'FontName','Calibri');

        if strcmp(vs,'prc')
            pctWet = 100*mean(y >= thrPrc); pctDry = 100 - pctWet;
            xline(ax, thrPrc, '--', 'LineWidth', 1.2);
            label = sprintf('Dry < %.1f mm: %.1f%%  |  Wet ≥ %.1f mm: %.1f%%', thrPrc, pctDry, thrPrc, pctWet);
        else
            medv = median(y,'omitnan'); pctLE = 100*mean(y <= medv); pctGT = 100 - pctLE;
            xline(ax, medv, '--', 'LineWidth', 1.2);
            label = sprintf('≤ Median: %.1f%%  |  > Median: %.1f%%', pctLE, pctGT);
        end
        text(ax, 0.02, 0.95, label, 'Units','normalized', 'HorizontalAlignment','left', ...
             'VerticalAlignment','top', 'FontWeight','bold', 'BackgroundColor',[1 1 1 0.75], 'Margin',4);

        ttl = sprintf('%s  (n = %d)', regs{r}, numel(y));
        title(ax, ttl, 'FontWeight','bold','Color', regColor.(regs{r}));
    end

    xlabel(tl, xLabelFA.(vs), 'FontName','Calibri','FontWeight','bold');
    title(tl, sprintf('منحنی چگالی احتمال ماهانه %s — %d–%d', upper(vs), Y0, Y1), 'FontWeight','bold');

    if ~exist(outDir,"dir"); mkdir(outDir); end
    outPNG = fullfile(outDir, sprintf('PDF_%s_%d_%d.png', vs, Y0, Y1));
    exportgraphics(tl, outPNG, 'Resolution', 200);
    fprintf('Saved: %s\n', outPNG);
end

function ncPath = findCruFile(dirCRU, varCRU)
    p1 = dir(fullfile(dirCRU, sprintf('cru_ts*.%s.dat.nc', varCRU)));
    if isempty(p1), p1 = dir(fullfile(dirCRU, sprintf('cru_ts*.%s.dat', varCRU))); end
    if isempty(p1), error('Could not find CRU file for %s in %s', varCRU, dirCRU); end
    [~,i] = max([p1.datenum]); ncPath = fullfile(p1(i).folder, p1(i).name);
end

function dt = cru_time_to_datetime(tnum, unitsStr)
    unitsStr = string(unitsStr);
    if contains(unitsStr, "months since")
        t0 = extractAfter(unitsStr, "months since "); t0 = datetime(char(t0),'InputFormat','yyyy-MM-dd');
        dt = dateshift(t0 + calmonths(tnum), 'start','month');
    elseif contains(unitsStr, "days since")
        t0 = extractAfter(unitsStr, "days since "); t0 = datetime(char(t0),'InputFormat','yyyy-MM-dd');
        dt = dateshift(t0 + days(tnum), 'start','month');
    else
        t0 = datetime(1900,1,1); dt = dateshift(t0 + calmonths(tnum), 'start','month');
    end
end

function [xi, fi] = compute_pdf(y)
    xi = []; fi = [];
    y = y(:); y = y(isfinite(y));
    if numel(y) < 5, return; end
    try
        lo = prctile(y,0.5); hi = prctile(y,99.5);
        if ~isfinite(lo) || ~isfinite(hi) || lo==hi, lo=min(y); hi=max(y); end
        xx = linspace(lo,hi,256);
        fi = ksdensity(y, xx, 'Function','pdf'); xi = xx;
    catch
        try
            [fi,edges] = histcounts(y,'Normalization','pdf');
            xi = 0.5*(edges(1:end-1)+edges(2:end));
        catch
            xi = []; fi = [];
        end
    end
end

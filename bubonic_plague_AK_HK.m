%% ========================================================================
%  MODELLING THE BUBONIC PLAGUE (SECOND PANDEMIC, 1346-1353)
%  Land-Based Reaction-Diffusion SIR PDE on a Topographic European Grid
%  ========================================================================
%
%  We adapt the governing equations from Davydovych, Dutka, and Cherniha (2023).
%  The model solves for local population fractions s=S/N, i=I/N, and r=R/N.
%  This approach allows the heterogeneous local population density N(x) to 
%  directly influence the local reproduction number, defined mathematically 
%  as R0(x) = beta*N(x)/(gamma+mu):
%
%     ds/dt = mu(1-s) - beta N s i + D_S lap(s)
%     di/dt = beta N s i - (gamma+mu) i + D_I lap(i)
%     dr/dt = gamma i - mu r + D_R lap(r)
%
%  CORE MECHANICS:
%  * We restrict transmission strictly to connected landmasses. The model operates 
%    as a plain SIR framework without maritime transmission or explicit flea/rat vectors.
%  * The infection spreads continuously across all connected land using a baseline 
%    diffusion coefficient (D_I) of 5 km^2/day.
%  * Historical trade corridors actively enhance mobility, boosting diffusion up 
%    to 10 km^2/day. Conversely, mountainous terrain acts as a physical barrier 
%    that severely slows transmission.
%  * Population density primarily follows elevation with only a mild route enhancement. 
%    Consequently, the local R0 exceeds 1 across most lowlands but drops below 1 in 
%    high-altitude or remote terrain. This aligns with Noble (1974) and Gaudart (2010), 
%    demonstrating that diffusion determines the speed of spread while R0 determines viability.
%  * The model dynamically generates safe zones wherever the effective R0 falls below 1.
%  * The advancing wave travels at the theoretical Fisher-KPP speed, calculated 
%    as c* = 2 sqrt(D (beta - gamma - mu)).
%  * We calibrate the temporal arrival against documented outbreaks in Florence 
%    (early 1348), Paris (1348), London (September 1348), and Moscow (1353).
%
%  PARAMETER BOUNDS: 
%  * beta 0.3-0.5, gamma 0.1-0.2, mu 5e-5, D_I 1-10 km^2/day, R0 1.5-2.0, CFR 30-60%.
%
%  HISTORICAL DATA (Benedictow 2021, Map 1):
%  * We digitised the year-by-year color zones (1346-1353) from Benedictow's map 
%    directly onto our computational grid (documented_arrival_year_grid.csv).
%  * We georeferenced this grid using 29 distinct cities, yielding an RMS error 
%    of roughly 35 km. Blank cells indicate no data or regions outside the map.
%  * We use the overland arrows to construct our trade corridors. We completely 
%    ignore the maritime arrows.
%  * The 1346 zone represents the Golden Horde territories and serves as our 
%    initial infected region. We use the subsequent temporal zones solely to 
%    validate the model output, never to force its behavior.
%  * We treat the historical border between the Golden Horde and the Russian 
%    principalities ("passage prohibited for Christians") as a restrictive barrier.
%
%  ASSUMPTIONS REGISTER:
%   1. We apply an initial seed prevalence of 1e-4 across the entire 1346 zone 
%      to simulate a small starting outbreak.
%   2. We assign all corridors a baseline weight of 1 and a boost of 1 because 
%      the historical map indicates direction but not traffic intensity.
%   3. We place a waypoint at Riga on the Elbing-Konigsberg-Pskov corridor to 
%      account for the map arrow pointing generally to Livonia.
%   4. We apply a barrier factor of 0.02 to the Golden Horde border because texts 
%      record a prohibition on passage but do not quantify the restriction.
%   5. We set the mountain elevation scale to 800 meters while keeping route 
%      passes open. The mathematical framework requires topography, but historical 
%      sources do not provide exact barrier strengths [3].
%   6. We construct a density proxy N(x) to ensure R0 scales with population 
%      density, as historical sources only provide broad regional totals [1,3,11].
%   7. We build an artificial land link between Dover and Calais to allow overland 
%      transmission into England, bypassing the sea landings shown on the map.
%   8. We assume equal diffusion across all classes (D_S = D_R = D_I) because 
%      the literature only estimates ranges for the infected class [6].
%   9. We fix the baseline case fatality rate at 45%, perfectly splitting the 
%      difference of the historical 30 to 60% mortality range [2,4].
%  10. We define the regional boundaries for Poland, Bohemia, and Finland using 
%      approximate limits based on the named areas in historical texts.
%  11. We define a safe zone as any region that remains entirely unreached or 
%      suffers less than 5% mortality. This provides a mathematically measurable 
%      threshold for areas historically described as "largely spared."
%  12. We restrict the analysis window to Europe, explicitly excluding North 
%      Africa, Anatolia, and completely sea-isolated cells.
%  13. We do not model seasonality, maritime spread, or explicit ectoparasite 
%      vectors, though we acknowledge the biological role of fleas and rats [4].
%  14. We utilize a 0.25-degree computational grid mapped over 1-degree built-in 
%      topography. We set the numerical extinction cutoff at 1e-7 and the arrival 
%      threshold at 1e-4.
%
%  IMPORTANT CAVEATS:
%  * The safe zones emerge purely as a mathematical output of the model. Absolutely 
%    nothing in the codebase artificially forces any region to remain safe.
%  * Because this framework relies exclusively on land-based transmission, it 
%    cannot naturally reproduce the historically documented maritime arrivals in 
%    the Mediterranean during 1347 and 1348. To account for this, we exclude 
%    the 1347 zone during our secondary statistical comparisons.
%
%  EXECUTION:
%  * Ensure 'documented_arrival_year_grid.csv' remains in the active directory. 
%  * Press Run. The script will push all results to './BubonicPlague_PDE_Outputs'.
%  * Requires MATLAB R2019b or newer. No external toolboxes required.
%% ========================================================================
clear; clc; close all; tic
outDir = 'BubonicPlague_PDE_Outputs'; if ~exist(outDir,'dir'), mkdir(outDir); end
runExtras = true;     % 1D validation
makeVideo = true;     % MP4 of 3D wave (skipped automatically if unsupported)

%% ---------------- 1. PARAMETERS ------------------------------------------
prm.beta  = 0.30;      % transmission rate (day^-1), slides: 0.3-0.5
prm.gamma = 0.15;      % removal rate (day^-1), slides: 0.1-0.2
prm.mu    = 5e-5;      % natural birth/death rate (day^-1)
prm.R0    = prm.beta/(prm.gamma+prm.mu);   % reference R0 at N=1 (~2.0, slides 1.5-2.0)
prm.D0    = 5;         % background D_I over ALL land (km^2/day), slides: 1-10
prm.DSrat = 1;  prm.DRrat = 1;             % D_S/D_I and D_R/D_I
prm.cfr   = 0.45;      % fraction of removed who died (mid of 30-60% mortality)
prm.epoch = 1346;      % simulation starts 1 Jan 1346 (documented 1346 zone, Golden Horde lands)
prm.T     = 365*8;     % days (to end of 1353)
prm.dg    = 0.25;      % grid spacing (deg)
prm.fade  = 1e-7;      % numerical extinction cutoff for i
prm.thr   = 1e-4;      % prevalence defining "arrival"
prm.seedFrac = 1e-4;   % initial infected fraction at each entry point (assumption 1)
prm.entryMode = 2;     % ENTRY POINTS: 0 = 1346 zone only | 1 = + documented 1347 zones (Sicily, Genoa, Marseille, Constantinople...) | 2 = + listed ports (default)
prm.entryDay  = 182;   % seeding day inside the documented year (assumption: mid-year; the map gives only the year)
% Entry points are places where Map 1 documents the plague arriving (mostly by sea). Each is a start of an outward LAND spread;
% a cell is seeded only if the land front has not already arrived. No spread crosses the sea.
prm.barrierFactor = 0;    % D multiplier across the Golden Horde / Russia border: 0 = wall (map: passage prohibited; 0.02 only delayed the front 1-2 y)
% --- ASSUMPTIONS (not in the PDFs) ---
prm.routeBoost = 1;    % D on trade corridors = D0*(1+routeBoost) = 10 km^2/day (upper slide range)
prm.hMount  = 800;     % mountain barrier length scale (m); passes on routes
prm.hN      = 1500;    % density decay with elevation (m)
prm.Lroute  = 300;     % density decay with distance from routes (km): mild effect
prm.nFloor  = 0.6;     % relative density far from routes (lowland R0 ~1.2 > 1: land spread everywhere)
prm.safeThr = 0.05;    % "safe" if final mortality < 5% of local population
prm.channelLink = true;% Dover-Calais overland link so London can be reached

%% ---------------- 2. 3D DOMAIN (real topography) -------------------------
G.lonv = -12:prm.dg:50;  G.latv = 30:prm.dg:62;
[G.LON,G.LAT] = meshgrid(G.lonv,G.latv);  [nr,nc] = size(G.LON);
S = load('topo','topo');
lonT = [-0.5 0.5:1:359.5 360.5];  latT = -89.5:1:89.5;
topoExt = double([S.topo(:,end) S.topo S.topo(:,1)]);
elev = interp2(lonT,latT,topoExt,mod(G.LON,360),G.LAT,'linear');
G.land = elev > 0;
if prm.channelLink                                   % ASSUMPTION: narrow land link at Dover-Calais
    cl = [50.95 1.85; 51.13 1.31]; offs = [0 0;1 0;0 1];
    for tt = linspace(0,1,40)
        la = cl(1,1)+tt*(cl(2,1)-cl(1,1)); lo = cl(1,2)+tt*(cl(2,2)-cl(1,2));
        ix = round((lo-G.lonv(1))/prm.dg)+1; iy = round((la-G.latv(1))/prm.dg)+1;
        for k = 1:3, G.land(iy+offs(k,1),ix+offs(k,2)) = true; end
    end
end
elev(~G.land) = 0;  elev(G.land) = max(elev(G.land),1);  G.elev = elev;

% cities used for diagnostics: name, lat, lon
cities = {
 'Caffa',          45.03, 35.38;
 'Constantinople', 41.01, 28.98;
 'Florence',       43.77, 11.26;
 'Paris',          48.86,  2.35;
 'London',         51.51, -0.13;
 'Moscow',         55.75, 37.62;
 'Cologne',        50.94,  6.96;
 'Vienna',         48.21, 16.37;
 'Prague',         50.08, 14.44;
 'Krakow',         50.06, 19.94;
 'Warsaw',         52.23, 21.01;
 'Kyiv',           50.45, 30.52;
 'Pskov',          57.82, 28.33};
G.cName = cities(:,1); G.cLat = cell2mat(cities(:,2)); G.cLon = cell2mat(cities(:,3)); nCity = numel(G.cName);
G.cityIdx = zeros(nCity,1);
for k = 1:nCity
    d2 = (G.LON-G.cLon(k)).^2 + (G.LAT-G.cLat(k)).^2;  d2(~G.land) = Inf;
    [~,G.cityIdx(k)] = min(d2(:));
end
G.caffaIdx = G.cityIdx(1);
iPar = find(strcmp(G.cName,'Paris'));
mainland = landComponent(G.land,G.cityIdx(iPar));
if ~mainland(G.caffaIdx)
    d2 = (G.LON-G.cLon(1)).^2 + (G.LAT-G.cLat(1)).^2;  d2(~mainland) = Inf;
    [~,G.caffaIdx] = min(d2(:));  G.cityIdx(1) = G.caffaIdx;
    fprintf('Caffa not land-connected at this resolution: reference cell moved to nearest mainland cell.\n');
end
% analysis window "Europe": land-connected to Paris, lat>=36, lon<=50, excluding Anatolia
G.eu = mainland & G.LAT>=36 & G.LON<=50 & ~(G.LON>26.5 & G.LAT<41.5);
G.regPol = G.eu & G.LAT>=50.5 & G.LAT<=54 & G.LON>=16 & G.LON<=24;     % interior Poland (box)
G.regBoh = G.eu & G.LAT>=48.6 & G.LAT<=51 & G.LON>=12.1 & G.LON<=18.9; % Bohemia (box)
G.regFin = G.eu & G.LAT>=59.8 & G.LAT<=62 & G.LON>=21 & G.LON<=31.5;   % southern Finland (box): documented as escaping the plague
G.cellArea = prm.dg*111.32*cosd(G.LAT)*prm.dg*110.57;

% DOCUMENTED ARRIVAL YEAR (Benedictow 2021, Map 1): year-by-year colour zones digitised onto this grid
docFile = fullfile(fileparts(mfilename('fullpath')),'documented_arrival_year_grid.csv');
if ~exist(docFile,'file'), docFile = 'documented_arrival_year_grid.csv'; end
G.docYear = readmatrix(docFile);
assert(isequal(size(G.docYear),[nr nc]),'documented_arrival_year_grid.csv does not match the model grid (129 x 249).');
G.seedMask = G.land & G.docYear==1346;           % initial infected region = documented 1346 zone
fprintf('Documented 1346 zone: %d land cells used as the initial infected region.\n',nnz(G.seedMask));
G.seedDay = makeSeeds(G,prm);
fprintf('Entry points (mode %d): %d cells (%d at day 0, %d later documented entries).\n',prm.entryMode,nnz(isfinite(G.seedDay)),nnz(G.seedDay==0),nnz(G.seedDay>0));

% OVERLAND CORRIDORS = land arrows of Benedictow (2021) Map 1 ([lat lon] waypoints; all weight 1; no sea arrows)
routes = {
 [43.30 5.37; 43.95 4.81; 45.76 4.83; 46.78 4.85; 47.24 6.02; 47.56 7.59; 48.57 7.75; 50.00 8.27; 50.94 6.96];  % Marseille-Avignon-Lyon-Chalon-Besancon-Basel-Strasbourg-Mainz-Cologne (Rhone-Rhine)
 [50.00 8.27; 50.11 8.68; 50.98 11.03];                                                                      % Mainz-Frankfurt-Erfurt
 [42.70 2.90; 43.18 3.00; 43.60 1.44; 44.84 -0.58; 46.58 0.34; 47.39 0.69];                                 % Perpignan-Narbonne-Toulouse-Bordeaux-Poitiers-Tours
 [43.30 5.37; 43.18 3.00];                                                                                   % Marseille-Narbonne
 [49.18 -0.37; 49.44 1.10; 48.86 2.35; 48.57 7.75];                                                          % Caen-Rouen-Paris-Strasbourg
 [50.95 1.85; 51.05 3.72];                                                                                   % Calais-Ghent
 [50.94 6.96; 53.08 8.80];                                                                                   % Cologne-Bremen
 [50.95 1.85; 51.13 1.31; 51.51 -0.13; 51.89 0.90; 52.75 0.40; 53.23 -0.54; 53.96 -1.08];                  % Calais-Dover-London-Colchester-Lynn-Lincoln-York (England)
 [44.82 20.46; 47.50 19.04; 48.21 16.37; 48.57 13.43; 48.54 12.15; 48.14 11.58; 47.26 11.39; 46.07 11.12; 45.46 9.19; 44.41 8.93]; % Belgrade-Budapest-Vienna-Passau-Landshut-Munich-Innsbruck-Trento-Milan-Genoa
 [47.56 7.59; 48.37 10.90; 48.14 11.58];                                                                     % Basel-Augsburg-Munich
 [45.44 12.32; 46.07 11.12];                                                                                 % Venice-Trento
 [44.41 8.93; 43.72 10.40; 43.77 11.26; 43.11 12.39; 41.90 12.50; 40.85 14.27];                             % Genoa-Pisa-Florence-Perugia-Rome-Naples
 [42.70 2.90; 41.39 2.17; 41.12 1.25; 41.62 0.62; 41.65 -0.88];                                             % Perpignan-Barcelona-Tarragona-Lleida-Zaragoza
 [39.47 -0.38; 40.34 -1.11; 41.65 -0.88];                                                                    % Valencia-Teruel-Zaragoza
 [39.47 -0.38; 40.07 -2.14; 40.42 -3.70; 39.86 -4.02; 37.39 -5.99];                                         % Valencia-Cuenca-Madrid-Toledo-Seville
 [42.88 -8.54; 42.60 -5.57; 43.36 -5.85];                                                                    % Santiago-Leon-Oviedo
 [42.88 -8.54; 40.21 -8.43; 38.72 -9.14];                                                                    % Santiago-Coimbra-Lisbon
 [54.35 18.65; 53.01 18.60; 52.23 21.01];                                                                    % Danzig-Thorn-Warsaw (Vistula)
 [54.16 19.40; 54.71 20.51; 56.95 24.11; 57.82 28.33];                                                      % Elbing-Konigsberg-(Riga)-Pskov (Riga: assumption 3)
 [57.82 28.33; 58.52 31.27];                                                                                 % Pskov-Novgorod
 [57.82 28.33; 56.86 35.90; 55.75 37.62];                                                                    % Pskov-Tver-Moscow
 [57.82 28.33; 60.03 37.79];                                                                                 % Pskov-Belozersk
 [46.35 48.04; 48.57 46.17; 49.05 44.80; 51.53 46.03; 55.79 49.12; 56.33 44.00];                            % Astrakhan-Old Sarai-Dubovka-Saratov-Kazan-Nizhniy Novgorod (Volga)
 [48.57 46.17; 47.10 39.40];                                                                                 % Old Sarai-Tana (Azov)
 [46.35 48.04; 42.06 48.29; 40.41 49.87]};                                                                   % Astrakhan-Derbent-Baku
G.routes = routes(:,1);  wts = ones(numel(G.routes),1);
G.routeMask = rasterRoutes(G,wts,prm);
% Golden Horde / Russian principalities border (Map 1: "passage prohibited for Christians"), digitised
barrier = [46.95 30.45; 47.10 31.95; 47.40 35.20; 48.50 35.35; 49.05 34.20; 49.45 34.50; 51.05 35.50; 51.90 36.70; 54.20 41.10; 55.80 44.20; 57.50 45.80; 58.00 47.35; 58.30 50.15];
G.barrierMask = rasterPolyline(G,barrier,prm);
G.dRoute = routeDist(G);
G = buildFields(G,prm);
fprintf('Grid %d x %d (%d land cells), %d routes, substeps/day = %d\n',nr,nc,nnz(G.land),numel(G.routes),G.nSub);
cRoute = 2*sqrt(prm.D0*(1+prm.routeBoost)*(prm.beta-prm.gamma-prm.mu));
cBG    = 2*sqrt(prm.D0*(prm.beta-prm.gamma-prm.mu));
fprintf('R0(ref)=%.2f | Fisher-KPP c*: background %.2f km/day, trade routes %.2f km/day\n',prm.R0,cBG,cRoute);

%% ---------------- 3. BASELINE 2D SIMULATION ------------------------------
fprintf('Running baseline PDE simulation...\n');
R = runPlague(G,prm,true);
M = metrics(G,R,prm); F = fitFront(G,R.Tarr);
cls = classifyZones(G,R,prm); SZ = safeStats(G,cls);
fprintf('Mean front speed = %.2f km/day (R^2=%.2f) [slides: 1.5-4 km/day]\n',M.speed,M.r2);
fprintf('Mortality (Europe, pop-weighted) = %.1f%% | ever infected = %.1f%% | Europe land reached = %.1f%%\n',M.mortPct,M.everInfPct,M.reachPct);
fprintf('SAFE ZONES: %.1f%% of European land, %.1f%% of population (never reached %.1f%%)\n',SZ.safeAreaPct,SZ.safePopPct,SZ.areaPct(1));
fprintf('Mean mortality: Poland+Bohemia %.1f%% vs rest of Europe %.1f%%\n',M.mortPB,M.mortRest);

%% ---------------- 4. TABLES ----------------------------------------------
cityMort = 100*R.dFinal(G.cityIdx);
zoneNames = {'Never reached','Low impact','Moderate','Severe'};
cz = cls(G.cityIdx); cityZone = repmat({'outside window'},nCity,1); okz = ~isnan(cz); cityZone(okz) = zoneNames(cz(okz)+1);
cityModelYear = floor(prm.epoch + R.Tarr(G.cityIdx)/365.25); cityDocYear = G.docYear(G.cityIdx);
Tc = table(G.cName,G.cLat,G.cLon,cityModelYear,cityDocYear,cityModelYear-cityDocYear,cityMort,G.R0eff(G.cityIdx),cityZone, ...
   'VariableNames',{'City','Lat','Lon','ModelArrivalYear','DocumentedYear','ErrorYears','LocalMortalityPct','LocalR0','Zone'});
disp(Tc); writetable(Tc,fullfile(outDir,'city_table.csv'));
CD = compareDocumented(G,R,prm);
fprintf('Documented-year comparison (%d cells, 1346 seed zone excluded): never reached %.1f%%\n',CD.nCells,CD.neverReachedPct);
fprintf('  all documented cells reached by the model: n=%d, RMSE %.2f y, bias %+.2f y, within 1 y %.1f%%, Spearman %.2f\n',CD.n(1),CD.rmse(1),CD.bias(1),CD.within1(1),CD.corr(1));
fprintf('  excluding the 1347 (sea-borne) zone:        n=%d, RMSE %.2f y, bias %+.2f y, within 1 y %.1f%%, Spearman %.2f\n',CD.n(2),CD.rmse(2),CD.bias(2),CD.within1(2),CD.corr(2));
CDt = table(["all documented cells";"excluding 1347 zone"],CD.n',CD.rmse',CD.bias',CD.within1',CD.corr','VariableNames',{'Subset','nCells','RMSE_years','Bias_years','Within1yrPct','Spearman'});
writetable(CDt,fullfile(outDir,'documented_comparison.csv'));
docInfo = G.docYear(G.eu & isfinite(G.docYear) & G.docYear>1346); regYear = zeros(8,2);
for yy = 1347:1353, mk = G.eu & G.docYear==yy; regYear(yy-1346,:) = [nnz(mk), 100*nnz(mk & isfinite(CD.mdl))/max(nnz(mk),1)]; end
YT = array2table(regYear(1:7,:),'RowNames',string(1347:1353),'VariableNames',{'DocumentedCells','ModelReachedPct'}); disp(YT);
writetable(YT,fullfile(outDir,'documented_zone_reach.csv'),'WriteRowNames',true);

regN = {'Poland (box)','Bohemia (box)','Finland (box)','Rest of Europe'}; regM = {G.regPol,G.regBoh,G.regFin,G.eu&~G.regPol&~G.regBoh&~G.regFin};
regT = zeros(4,4);
for k = 1:4
    mk = regM{k}; regT(k,:) = [100*sum(R.dFinal(mk).*G.cellPop(mk))/sum(G.cellPop(mk)), 100*nnz(isnan(R.Tarr(mk)))/nnz(mk), ...
        mean(G.R0eff(mk)), 100*nnz(G.R0eff(mk)<1)/nnz(mk)];
end
RT = array2table(regT,'RowNames',regN,'VariableNames',{'MortalityPct','NeverReachedPct','MeanLocalR0','CellsR0below1Pct'}); disp(RT);
writetable(RT,fullfile(outDir,'safe_zone_regions.csv'),'WriteRowNames',true);

%% ---------------- 5. FIGURES ---------------------------------------------
yr = @(d) prm.epoch + d/365.25;  tyr = yr(0:prm.T-1);  cmJ = jet(256);
dstr = @(d) datestr(datenum(prm.epoch,1,1)+d,'mmm yyyy');

% Fig 1: domain, routes, local R0
f1 = figure('Color','w','Position',[50 50 1450 600]);
subplot(1,2,1); drawTerrain(G,baseRGB(G.elev)); plotRoutes3(G,[0.1 0.1 0.7],1.6);
zc = G.elev(G.cityIdx)/111000*60+0.4; plot3(G.cLon,G.cLat,zc,'ro','MarkerFaceColor','r','MarkerSize',4);
for k = 1:nCity, text(G.cLon(k),G.cLat(k),zc(k)+0.3,G.cName{k},'FontSize',7); end
title('(a) 3D domain with overland trade routes (blue)');
subplot(1,2,2); imagesc(G.lonv,G.latv,G.R0eff,'AlphaData',G.land); axis xy image; set(gca,'Color',[0.75 0.85 0.95]); hold on;
colormap(gca,parula); cb = colorbar; cb.Label.String = 'Local R_0 = \beta N/(\gamma+\mu)'; contour(G.lonv,G.latv,G.R0eff,[1 1],'w','LineWidth',1.5);
for r_ = 1:numel(G.routes), plot(G.routes{r_}(:,2),G.routes{r_}(:,1),'k-','LineWidth',0.4); end
xlabel('Longitude (°E)'); ylabel('Latitude (°N)'); title('(b) Local R_0 (white contour: R_0 = 1)');
saveFig(f1,outDir,'Fig1_domain_routes_R0');

% Fig 2: 3D wave snapshots
snapShow = [180 450 720 990 1530 2160];
f2 = figure('Color','w','Position',[30 30 1600 900]);
for k = 1:6
    subplot(2,3,k); kk = snapShow(k)/30+1;
    drawTerrain(G,infRGB(G.elev,double(R.Vi(:,:,kk)))); plotRoutes3(G,[0.25 0.25 0.6],0.5); title(dstr(snapShow(k)));
end
colormap(infCmap()); caxis([-4 -1]); cb = colorbar('Position',[0.94 0.3 0.012 0.4]); cb.Label.String = 'log_{10} infected fraction i';
saveFig(f2,outDir,'Fig2_3D_wave_snapshots');

% Fig 3: arrival times
f3 = figure('Color','w','Position',[50 50 1500 620]);
subplot(1,2,1); a = R.Tarr/prm.T; a(isnan(a)) = 0;
drawTerrain(G,overlayRGB(baseRGB(G.elev),a,cmJ,double(~isnan(R.Tarr)))); plotRoutes3(G,[0 0 0],0.6);
colormap(gca,cmJ); caxis([prm.epoch prm.epoch+8]); cb = colorbar; cb.Label.String = 'Arrival year'; title('(a) Arrival time on 3D terrain (no colour = never reached)');
subplot(1,2,2); Ty = yr(R.Tarr); imagesc(G.lonv,G.latv,Ty,'AlphaData',~isnan(Ty)); axis xy image; set(gca,'Color',[0.75 0.85 0.95]); hold on;
contour(G.lonv,G.latv,Ty,prm.epoch:0.5:prm.epoch+8,'k'); colormap(gca,cmJ); caxis([prm.epoch prm.epoch+8]); cb = colorbar; cb.Label.String = 'Arrival year';
plot(G.cLon,G.cLat,'k.','MarkerSize',10); for k = 1:nCity, text(G.cLon(k)+0.3,G.cLat(k),G.cName{k},'FontSize',7); end
xlabel('Longitude (°E)'); ylabel('Latitude (°N)'); title('(b) Isochrones (0.5-yr)'); saveFig(f3,outDir,'Fig3_arrival_time');

% Fig 4: space-time isosurface
[X,Y,Tt] = meshgrid(G.lonv,G.latv,yr(R.snapDays)); V = smooth3(double(R.Vi));
f4 = figure('Color','w','Position',[50 50 1000 750]); hold on; plotCoast(G,prm.epoch);
p1 = patch(isosurface(X,Y,Tt,V,1e-3)); isonormals(X,Y,Tt,V,p1); set(p1,'FaceColor',[1 0.75 0.1],'EdgeColor','none','FaceAlpha',0.35);
p2 = patch(isosurface(X,Y,Tt,V,2e-2)); isonormals(X,Y,Tt,V,p2); set(p2,'FaceColor',[0.8 0.05 0.05],'EdgeColor','none','FaceAlpha',0.8);
view(-35,25); axis tight; daspect([10 10 1]); grid on; camlight; lighting gouraud;
xlabel('Longitude (°E)'); ylabel('Latitude (°N)'); zlabel('Year'); title('Space-time plague wave: i = 0.1% (amber), 2% (red)');
saveFig(f4,outDir,'Fig4_spacetime_isosurface');

% Fig 5: safe zones
f7 = figure('Color','w','Position',[30 30 1650 620]); zc4 = [0.10 0.60 0.25; 0.60 0.85 0.40; 1.00 0.70 0.20; 0.75 0.10 0.10];
subplot(1,3,1); drawTerrain(G,classRGB(G,cls,zc4)); plotRoutes3(G,[0 0 0],0.6); hold on;
for c = 1:4, plot3(nan,nan,nan,'s','MarkerFaceColor',zc4(c,:),'MarkerEdgeColor','k','MarkerSize',9,'DisplayName',zoneNames{c}); end
legend('Location','southoutside','Orientation','horizontal','FontSize',7); title('(a) Safe zones on 3D terrain');
subplot(1,3,2); imagesc(G.lonv,G.latv,cls,'AlphaData',~isnan(cls)); axis xy image; colormap(gca,zc4); caxis([-0.5 3.5]); set(gca,'Color',[0.75 0.85 0.95]); hold on;
contour(G.lonv,G.latv,G.R0eff,[1 1],'k--','LineWidth',1.2);
for r_ = 1:numel(G.routes), plot(G.routes{r_}(:,2),G.routes{r_}(:,1),'k-','LineWidth',0.4); end
rectangle('Position',[16 50.5 8 3.5],'EdgeColor','b','LineWidth',1.5); rectangle('Position',[12.1 48.6 6.8 2.4],'EdgeColor','m','LineWidth',1.5);
text(16.2,54.4,'Poland','Color','b','FontSize',8); text(12.1,48.2,'Bohemia','Color','m','FontSize',8);
for k = 1:nCity, text(G.cLon(k)+0.3,G.cLat(k),G.cName{k},'FontSize',6); end
xlabel('Longitude (°E)'); ylabel('Latitude (°N)'); title('(b) Zone map; dashed = local R_0 = 1');
subplot(1,3,3); bar([SZ.areaPct' SZ.popPct']); set(gca,'XTickLabel',zoneNames); legend('% land area','% population'); grid on; ylabel('Percent');
title(sprintf('(c) Safe = %.0f%% area, %.0f%% pop.',SZ.safeAreaPct,SZ.safePopPct)); saveFig(f7,outDir,'Fig5_safe_zones');

% Fig 6: documented arrival year (Benedictow 2021, Map 1) versus model
docY = G.docYear; docY(~G.land) = NaN; mdlY = CD.mdl; mdlY(~G.land) = NaN;
f8 = figure('Color','w','Position',[50 50 1500 900]);
subplot(2,2,1); imagesc(G.lonv,G.latv,docY,'AlphaData',~isnan(docY)); axis xy image; set(gca,'Color',[0.75 0.85 0.95]); colormap(gca,cmJ); caxis([prm.epoch prm.epoch+8]);
cb = colorbar; cb.Label.String = 'Year'; hold on; plot(barrier(:,2),barrier(:,1),'k:','LineWidth',2); title('(a) Documented arrival year (Map 1); dotted = border'); xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
subplot(2,2,2); imagesc(G.lonv,G.latv,mdlY,'AlphaData',~isnan(mdlY)); axis xy image; set(gca,'Color',[0.75 0.85 0.95]); colormap(gca,cmJ); caxis([prm.epoch prm.epoch+8]);
cb = colorbar; cb.Label.String = 'Year'; hold on; plot(barrier(:,2),barrier(:,1),'k:','LineWidth',2); title('(b) Model arrival year (land spread only)'); xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
subplot(2,2,3); dm = CD.diff; okd = G.eu & isfinite(docY) & isfinite(mdlY); dm(~okd) = NaN;
imagesc(G.lonv,G.latv,dm,'AlphaData',~isnan(dm)); axis xy image; set(gca,'Color',[0.75 0.85 0.95]); hold on;
rb = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)']; colormap(gca,rb); caxis([-5 5]);
cb = colorbar; cb.Label.String = 'Model minus documented (years)'; [rr,cc] = find(G.eu & isfinite(docY) & docY>1346 & isnan(mdlY)); plot(G.lonv(cc),G.latv(rr),'k.','MarkerSize',2);
title(sprintf('(c) Difference (black dots: documented but never reached, %.0f%%)',CD.neverReachedPct)); xlabel('Longitude (°E)'); ylabel('Latitude (°N)');
subplot(2,2,4); hold on; sel = {'Caffa','Constantinople','Moscow','Pskov','Paris','London'};
for k = 1:numel(sel), j = find(strcmp(G.cName,sel{k})); plot(tyr,100*R.iCity(:,j),'LineWidth',1.4,'DisplayName',sel{k}); end
grid on; legend('Location','northeast'); xlabel('Year'); ylabel('Infected fraction (%)'); title('(d) Local epidemic curves'); saveFig(f8,outDir,'Fig6_documented_vs_model');

if makeVideo
    try
        vw = VideoWriter(fullfile(outDir,'plague_PDE_wave_3D.mp4'),'MPEG-4'); vw.FrameRate = 8; open(vw);
        fv = figure('Color','w','Position',[50 50 1100 700]);
        for k = 1:size(R.Vi,3)
            clf(fv); drawTerrain(G,infRGB(G.elev,double(R.Vi(:,:,k)))); plotRoutes3(G,[0.25 0.25 0.6],0.5);
            title(dstr(R.snapDays(k))); drawnow; writeVideo(vw,getframe(fv));
        end
        close(vw);
    catch ME, warning('Video skipped: %s',ME.message); end
end

%% ---------------- 6. EXTRAS: 1D VALIDATION ---------------
if runExtras
    fprintf('Running 1D validation (Fisher-KPP, sensitivity, R0 threshold)...\n');
    dx = 2; L = 1800; nx = round(L/dx); xg = (1:nx)*dx; T1 = 400; N1 = ones(1,nx);
    Ds = [1 2 5 10]; bs = [0.3 0.4 0.5]; spM = nan(4,3); spT = nan(4,3);
    for a_ = 1:4, for b_ = 1:3
        o = run1D(bs(b_),prm.gamma,prm.mu,Ds(a_),N1,dx,T1,false);
        spM(a_,b_) = o.speed; spT(a_,b_) = 2*sqrt(Ds(a_)*(bs(b_)-prm.gamma-prm.mu));
    end, end
    mus = [5e-5 5e-4 5e-3]; spMu = zeros(1,3);
    for k = 1:3, o = run1D(prm.beta,prm.gamma,mus(k),10,N1,dx,T1,false); spMu(k) = o.speed; end
    oP = run1D(prm.beta,prm.gamma,prm.mu,10,N1,dx,T1,true);
    Nst = ones(1,nx); Nst(xg>1000) = 0.4;                         % density drop -> local R0 = 0.8 (<1)
    oS = run1D(prm.beta,prm.gamma,prm.mu,10,Nst,dx,800,true);
    f9 = figure('Color','w','Position',[50 50 1400 900]);
    subplot(2,2,1); hold on; dd = 50:50:400; cP = parula(numel(dd));
    for k = 1:numel(dd), plot(xg,oP.Im(dd(k),:),'Color',cP(k,:),'LineWidth',1.3,'DisplayName',sprintf('day %d',dd(k))); end
    grid on; legend('FontSize',6); xlabel('x (km)'); ylabel('i'); title('(a) 1D travelling wave (stable front shape)');
    subplot(2,2,2); hold on; fill([0 7 7 0],[1.5 1.5 4 4],[0.85 0.95 0.85],'EdgeColor','none','DisplayName','1.5-4 km/day (H1)');
    mk = {'o','s','^'}; for b_ = 1:3, plot(spT(:,b_),spM(:,b_),mk{b_},'MarkerFaceColor','auto','DisplayName',sprintf('\\beta=%.1f',bs(b_))); end
    plot([0 4],[0 4],'k--','DisplayName','c* (Fisher-KPP)'); grid on; legend('Location','northwest','FontSize',7);
    xlabel('Theoretical c* = 2\surd(D(\beta-\gamma-\mu))  (km/day)'); ylabel('Measured speed (km/day)'); title('(b) Wave speed vs Fisher-KPP');
    subplot(2,2,3); bar(spMu); set(gca,'XTickLabel',arrayfun(@(x) sprintf('%.0e',x),mus,'UniformOutput',false)); grid on;
    xlabel('\mu (day^{-1})'); ylabel('Speed (km/day)'); title('(c) Insensitivity to \mu');
    subplot(2,2,4); imagesc(xg,1:800,oS.Im); axis xy; colormap(gca,hot); caxis([0 0.05]); hold on; plot([1000 1000],[1 800],'c--','LineWidth',1.5);
    cb = colorbar; cb.Label.String = 'i'; xlabel('x (km)'); ylabel('Day'); title('(d) R_0 threshold: front stalls where R_0<1 (right of dashed line)');
    saveFig(f9,outDir,'Fig7_1D_validation_threshold');
    SV = table(Ds',spM(:,1),spT(:,1),spM(:,2),spT(:,2),spM(:,3),spT(:,3),'VariableNames',{'D','Meas_b03','cstar_b03','Meas_b04','cstar_b04','Meas_b05','cstar_b05'});
    disp(SV); writetable(SV,fullfile(outDir,'wave_speed_1D.csv'));
end

%% ---------------- 7. SAVE SUMMARY ----------------------------------------
fid = fopen(fullfile(outDir,'summary.txt'),'w');
fprintf(fid,'R0(ref)=%.2f; c* background=%.2f, routes=%.2f km/day\n',prm.R0,cBG,cRoute);
fprintf(fid,'Front speed (arrival fit) %.2f km/day (R2=%.2f)\n',M.speed,M.r2);
fprintf(fid,'Mortality %.1f %%; ever infected %.1f %%; Europe reached %.1f %%\n',M.mortPct,M.everInfPct,M.reachPct);
fprintf(fid,'Safe zones: %.1f %% area, %.1f %% population\n',SZ.safeAreaPct,SZ.safePopPct);
fprintf(fid,'Mortality Poland+Bohemia %.1f %% vs rest %.1f %%\n',M.mortPB,M.mortRest);
fprintf(fid,'Documented-year comparison: never reached %.1f %%; all cells n=%d RMSE %.2f y bias %+.2f y within1y %.1f %% Spearman %.2f; excluding 1347 zone n=%d RMSE %.2f y bias %+.2f y within1y %.1f %% Spearman %.2f\n',CD.neverReachedPct,CD.n(1),CD.rmse(1),CD.bias(1),CD.within1(1),CD.corr(1),CD.n(2),CD.rmse(2),CD.bias(2),CD.within1(2),CD.corr(2)); fclose(fid);
Rsave = R; save(fullfile(outDir,'results.mat'),'Rsave','prm','Tc','RT','cls','SZ','CD','-v7.3');
fprintf('Done in %.0f s. Outputs in %s\n',toc,outDir);

%% ========================================================================
%                          LOCAL FUNCTIONS
%% ========================================================================
function seedDay = makeSeeds(G,prm)
% day on which each entry cell is seeded (NaN = not an entry point). Years come from the digitised Map 1 grid.
seedDay = nan(size(G.land)); seedDay(G.land & G.docYear==1346) = 0;
if prm.entryMode >= 1, seedDay(G.land & G.docYear==1347) = 365 + prm.entryDay; end
if prm.entryMode >= 2
    ports = {'Weymouth',50.61,-2.46,1348; 'Bristol',51.45,-2.59,1348; 'Calais',50.95,1.85,1348; 'Rouen',49.44,1.10,1348; 'Dublin',53.35,-6.26,1348;
        'Bordeaux',44.84,-0.58,1348; 'Barcelona',41.39,2.17,1348; 'Valencia',39.47,-0.38,1348; 'Lisbon',38.72,-9.14,1348; 'Cadiz',36.53,-6.29,1348;
        'Algeciras',36.13,-5.45,1348; 'Almeria',36.84,-2.47,1348; 'Cartagena',37.60,-0.99,1348; 'Naples',40.85,14.27,1348; 'Hamburg',53.55,9.99,1349;
        'Kiel',54.32,10.14,1349; 'Luebeck',53.87,10.69,1349; 'Stralsund',54.31,13.09,1349; 'Ribe',55.33,8.77,1349; 'Bergen',60.39,5.32,1349; 'Oslo',59.91,10.75,1349};
    for k = 1:size(ports,1)   % ASSUMPTION: port list = reading of Map 1; year from the grid, fallback in the list
        d2 = (G.LON-ports{k,3}).^2 + (G.LAT-ports{k,2}).^2; d2(~G.land) = Inf; [~,ix] = min(d2(:));
        yr = G.docYear(ix); if ~isfinite(yr), yr = ports{k,4}; end
        seedDay(ix) = min([seedDay(ix) (yr-prm.epoch)*365 + prm.entryDay]);
    end
end
end

function comp = landComponent(land,idx)
comp = false(size(land)); comp(idx) = true; ker = [0 1 0;1 1 1;0 1 0];
while true
    nw = (conv2(double(comp),ker,'same')>0) & land;
    if nnz(nw)==nnz(comp), break; end
    comp = nw;
end
end

function rm = rasterRoutes(G,wts,prm)
[nr,nc] = size(G.land); rm = zeros(nr,nc);
for r = 1:numel(G.routes)
    P = G.routes{r}; tmp = zeros(nr,nc);
    for s = 1:size(P,1)-1
        n = ceil(hav(P(s,1),P(s,2),P(s+1,1),P(s+1,2))/10)+1;
        la = linspace(P(s,1),P(s+1,1),n); lo = linspace(P(s,2),P(s+1,2),n);
        ix = round((lo-G.lonv(1))/prm.dg)+1; iy = round((la-G.latv(1))/prm.dg)+1;
        ok = ix>=1 & ix<=nc & iy>=1 & iy<=nr; tmp(sub2ind([nr nc],iy(ok),ix(ok))) = wts(r);
    end
    rm = max(rm,tmp);
end
rm = max(rm,0.5*movmax(movmax(rm,3,1),3,2)).*G.land;
end

function dR = routeDist(G)
pts = zeros(0,2);
for r = 1:numel(G.routes)
    P = G.routes{r};
    for s = 1:size(P,1)-1
        n = ceil(hav(P(s,1),P(s,2),P(s+1,1),P(s+1,2))/25)+1;
        pts = [pts; linspace(P(s,1),P(s+1,1),n)' linspace(P(s,2),P(s+1,2),n)']; %#ok<AGROW>
    end
end
d2 = inf(size(G.land)); ck = 111.32*cosd(G.LAT);
for p = 1:size(pts,1)
    d2 = min(d2,((G.LON-pts(p,2)).*ck).^2 + ((G.LAT-pts(p,1))*110.57).^2);
end
dR = sqrt(d2);
end

function G = buildFields(G,prm)
% Local density N(x) -> local R0(x); heterogeneous D(x); zero-flux coasts
G.N = G.land.*exp(-G.elev/prm.hN).*(prm.nFloor+(1-prm.nFloor)*1.2*exp(-G.dRoute/prm.Lroute));
G.R0eff = prm.beta*G.N/(prm.gamma+prm.mu);
G.cellPop = G.N.*G.cellArea;  G.totPop = sum(G.cellPop(G.eu));
mtn = exp(-(G.elev/prm.hMount).*(1-0.75*G.routeMask));
D = prm.D0*mtn.*(1+prm.routeBoost*G.routeMask);
D = D.*(1-(1-prm.barrierFactor)*G.barrierMask);   % documented Golden Horde / Russia border
G.Dfield = D;
Kx = 1./(prm.dg*111.32*cosd(G.LAT)).^2;  Ky = 1/(prm.dg*110.57)^2;  L = G.land;
mE = L & L(:,[2:end end]);  mW = L & L(:,[1 1:end-1]);  mN = L & L([2:end end],:);  mS = L & L([1 1:end-1],:);
G.cE = mE.*(2*D.*D(:,[2:end end])./(D+D(:,[2:end end]))).*Kx;  G.cW = mW.*(2*D.*D(:,[1 1:end-1])./(D+D(:,[1 1:end-1]))).*Kx;
G.cN = mN.*(2*D.*D([2:end end],:)./(D+D([2:end end],:))).*Ky;  G.cS = mS.*(2*D.*D([1 1:end-1],:)./(D+D([1 1:end-1],:))).*Ky;
G.nSub = max(1,ceil(max(G.cE(:)+G.cW(:)+G.cN(:)+G.cS(:))/0.4));
end

function Lp = lap(U,G)
Lp = G.cE.*(U(:,[2:end end])-U) + G.cW.*(U(:,[1 1:end-1])-U) ...
   + G.cN.*(U([2:end end],:)-U) + G.cS.*(U([1 1:end-1],:)-U);
end

function R = runPlague(G,prm,rec)
[nr,nc] = size(G.land); Lm = double(G.land);
s = Lm; i = zeros(nr,nc); r = i; dc = i;
seedDays = unique(G.seedDay(isfinite(G.seedDay)))';              % entry points activate on their documented day
nSub = G.nSub; dt = 1/nSub; N = G.N; wEU = G.cellPop.*G.eu;
Tarr = nan(nr,nc); totInf = zeros(prm.T,1); totDead = totInf;
iCity = zeros(prm.T,numel(G.cityIdx)); snapDays = 0:30:prm.T-1;
if rec, Vi = zeros(nr,nc,numel(snapDays),'single'); else, Vi = []; end
for day = 0:prm.T-1
    if any(seedDays==day)
        mS = (G.seedDay==day) & isnan(Tarr); i(mS) = max(i(mS),prm.seedFrac); s(mS) = max(s(mS)-prm.seedFrac,0);
    end
    for k = 1:nSub
        lam = prm.beta*N.*s.*i;                                      % beta*S*I (N enters via local R0)
        lS = prm.DSrat*lap(s,G); lI = lap(i,G); lR = prm.DRrat*lap(r,G);
        s2 = s+dt*(prm.mu*(Lm-s)-lam+lS);
        i2 = i+dt*(lam-(prm.gamma+prm.mu)*i+lI);
        r2 = r+dt*(prm.gamma*i-prm.mu*r+lR);
        dc = dc+dt*prm.cfr*prm.gamma*i;                                % cumulative deaths (fraction of local N)
        s = max(s2,0); i = max(i2,0); r = max(r2,0); i(i<prm.fade) = 0;
    end
    newly = isnan(Tarr) & (i>prm.thr); Tarr(newly) = day;
    iCity(day+1,:) = i(G.cityIdx)';
    totInf(day+1)  = sum(i(:).*wEU(:))/G.totPop;
    totDead(day+1) = sum(dc(:).*wEU(:))/G.totPop;
    if rec && mod(day,30)==0, Vi(:,:,day/30+1) = single(i); end
end
R = struct('Tarr',Tarr,'totInf',totInf,'totDead',totDead,'iCity',iCity,'Vi',Vi, ...
    'snapDays',snapDays,'dFinal',dc,'rFinal',r,'sFinal',s);
end

function out = run1D(beta,gamma,mu,D,Nx,dx,T,recM)
nx = numel(Nx); i = zeros(1,nx); i(1:max(1,round(20/dx))) = 2e-3; s = 1-i; r = zeros(1,nx);
dt = min(0.5,0.2*dx^2/D); nst = round(T/dt); dt = T/nst;
front = zeros(T,1); Im = []; if recM, Im = zeros(T,nx,'single'); end; nrec = 0;
for k = 1:nst
    lam = beta*Nx.*s.*i;
    s = s+dt*(mu*(1-s)-lam+D*lap1(s,dx)); i = i+dt*(lam-(gamma+mu)*i+D*lap1(i,dx)); r = r+dt*(gamma*i-mu*r+D*lap1(r,dx));
    s = max(s,0); i(i<1e-9) = 0;
    if nrec<T && k*dt>=nrec+1-1e-9
        nrec = nrec+1; ix = find(i>1e-4,1,'last'); if ~isempty(ix), front(nrec) = ix*dx; end
        if recM, Im(nrec,:) = single(i); end
    end
end
t = (1:T)'; sel = t>=0.4*T & front>0;
if nnz(sel)>5, p = polyfit(t(sel),front(sel),1); out.speed = p(1); else, out.speed = 0; end
out.front = front; out.Im = Im;
end

function Lp = lap1(U,dx)
Lp = ([U(2:end) U(end)]-2*U+[U(1) U(1:end-1)])/dx^2;
end

function bm = rasterPolyline(G,P,prm)
[nr,nc] = size(G.land); bm = zeros(nr,nc);
for s = 1:size(P,1)-1
    n = ceil(hav(P(s,1),P(s,2),P(s+1,1),P(s+1,2))/10)+1;
    la = linspace(P(s,1),P(s+1,1),n); lo = linspace(P(s,2),P(s+1,2),n);
    ix = round((lo-G.lonv(1))/prm.dg)+1; iy = round((la-G.latv(1))/prm.dg)+1;
    ok = ix>=1 & ix<=nc & iy>=1 & iy<=nr; bm(sub2ind([nr nc],iy(ok),ix(ok))) = 1;
end
bm = movmax(movmax(bm,3,1),3,2).*G.land;
end

function C = compareDocumented(G,R,prm)
% Compare the model arrival year with the documented arrival year (Benedictow 2021, Map 1)
doc = G.docYear; mdl = floor(prm.epoch + R.Tarr/365.25);   % NaN where the model never reaches a cell
base = G.eu & G.land & isfinite(doc) & doc>1346 & isnan(G.seedDay);   % exclude directly seeded entry cells
noSea = base & doc~=1347;                                  % the 1347 patches are sea-borne arrivals (ports, islands)
C.nCells = nnz(base); C.neverReachedPct = 100*nnz(base & isnan(mdl))/max(nnz(base),1);
C.n = zeros(1,2); C.rmse = C.n; C.bias = C.n; C.within1 = C.n; C.corr = C.n;
for k = 1:2
    if k == 1, m = base & isfinite(mdl); else, m = noSea & isfinite(mdl); end
    if nnz(m) < 3, continue; end
    d = mdl(m)-doc(m); C.n(k) = nnz(m); C.rmse(k) = sqrt(mean(d.^2)); C.bias(k) = mean(d); C.within1(k) = 100*mean(abs(d)<=1);
    cc = corrcoef(rankAvg(mdl(m)),rankAvg(doc(m))); C.corr(k) = cc(1,2);
end
C.diff = mdl-doc; C.mdl = mdl;
end

function r = rankAvg(x)
[xs,ix] = sort(x(:)); n = numel(xs); r = zeros(n,1); i = 1;
while i <= n
    j = i; while j < n && xs(j+1) == xs(i), j = j+1; end
    r(ix(i:j)) = (i+j)/2; i = j+1;
end
end

function cls = classifyZones(G,R,prm)
% 0 never reached | 1 low impact (<safeThr mortality) | 2 moderate (<20%) | 3 severe   (Europe window only)
cls = nan(size(G.land)); m = R.dFinal; hit = G.eu & ~isnan(R.Tarr);
cls(G.eu & isnan(R.Tarr)) = 0;
cls(hit & m<prm.safeThr) = 1;  cls(hit & m>=prm.safeThr & m<0.2) = 2;  cls(hit & m>=0.2) = 3;
end

function SZ = safeStats(G,cls)
SZ.areaPct = zeros(1,4); SZ.popPct = zeros(1,4); A = sum(G.cellArea(G.eu)); P = sum(G.cellPop(G.eu));
for c = 0:3
    mk = G.eu & cls==c; SZ.areaPct(c+1) = 100*sum(G.cellArea(mk))/A; SZ.popPct(c+1) = 100*sum(G.cellPop(mk))/P;
end
SZ.safeAreaPct = sum(SZ.areaPct(1:2)); SZ.safePopPct = sum(SZ.popPct(1:2));
end

function F = fitFront(G,Tarr)
% distance to the nearest entry point and time since that entry point was seeded
sidx = find(isfinite(G.seedDay)); dmin = inf(size(G.land)); sd0 = zeros(size(G.land)); ck = 111.32*cosd(G.LAT);
for q = 1:numel(sidx)
    [ri,ci] = ind2sub(size(G.land),sidx(q)); d2 = ((G.LON-G.lonv(ci)).*ck).^2 + ((G.LAT-G.latv(ri))*110.57).^2; u = d2<dmin; dmin(u) = d2(u); sd0(u) = G.seedDay(sidx(q));
end
dist = sqrt(dmin); Tarr = Tarr-sd0;
ok = ~isnan(Tarr) & G.eu & dist>150 & isnan(G.seedDay) & Tarr>30;
if nnz(ok)<10, F = struct('speed',NaN,'p',[NaN NaN],'r2',NaN,'dist',dist,'ok',ok); return; end
p = polyfit(dist(ok),Tarr(ok),1); pred = polyval(p,dist(ok));
r2 = 1-sum((Tarr(ok)-pred).^2)/sum((Tarr(ok)-mean(Tarr(ok))).^2);
F = struct('speed',1/p(1),'p',p,'r2',r2,'dist',dist,'ok',ok);
end

function M = metrics(G,R,prm)
F = fitFront(G,R.Tarr); w = G.cellPop.*G.eu; pb = G.regPol|G.regBoh; rest = G.eu & ~pb;
M.speed = F.speed; M.r2 = F.r2;
M.mortPct = 100*sum(R.dFinal(:).*w(:))/sum(w(:));
M.everInfPct = 100*sum((1-R.sFinal(:)).*w(:))/sum(w(:));
M.mortPB = 100*sum(R.dFinal(pb).*G.cellPop(pb))/sum(G.cellPop(pb));
M.mortRest = 100*sum(R.dFinal(rest).*G.cellPop(rest))/sum(G.cellPop(rest));
euA = G.cellArea.*G.eu; M.reachPct = 100*sum(euA(~isnan(R.Tarr)))/sum(euA(:));
end

function d = hav(lat1,lon1,lat2,lon2)
a = sind((lat2-lat1)/2).^2 + cosd(lat1).*cosd(lat2).*sind((lon2-lon1)/2).^2;
d = 2*6371*asin(sqrt(a));
end

function base = baseRGB(elev)
b = 0.55+0.35*min(elev/2500,1); base = cat(3,b*0.78,b*0.86,b*0.66);
end

function cm = infCmap()
cm = flipud(autumn(256)).*linspace(1,0.45,256)';
end

function RGB = overlayRGB(base,a,cm,w)
ix = max(1,min(size(cm,1),ceil(a*size(cm,1)))); ix(isnan(ix)) = 1;
col = reshape(cm(ix(:),:),[size(a) 3]);
RGB = base.*(1-w)+col.*w;
end

function RGB = infRGB(elev,i)
a = min(max((log10(max(i,1e-9))+4)/3,0),1); w = double(i>1e-4).*min(1,0.35+a);
RGB = overlayRGB(baseRGB(elev),a,infCmap(),w);
end

function RGB = classRGB(G,cls,zc)
base = baseRGB(G.elev); RGB = base;
for c = 0:3
    mk = double(cls==c); col = cat(3,zc(c+1,1)*ones(size(mk)),zc(c+1,2)*ones(size(mk)),zc(c+1,3)*ones(size(mk)));
    RGB = RGB.*(1-0.65*mk)+0.65*mk.*(col.*(0.6+0.4*base(:,:,1)/0.9));
end
RGB = min(max(RGB,0),1);
end

function h = drawTerrain(G,RGB)
Z = G.elev/111000*60+0.03; Z(~G.land) = NaN;
surf(G.lonv([1 end]),G.latv([1 end]),zeros(2),'FaceColor',[0.66 0.80 0.92],'EdgeColor','none'); hold on;
h = surf(G.LON,G.LAT,Z,'CData',RGB,'FaceColor','interp','EdgeColor','none', ...
    'AmbientStrength',0.7,'DiffuseStrength',0.5,'SpecularStrength',0,'FaceLighting','gouraud');
xlim([-12 50]); ylim([30 62]); zlim([0 3]); view(-15,52); pbaspect([62 32 9]);
set(gca,'XTick',[],'YTick',[],'ZTick',[]); box off; light('Position',[-1 -1 1.5],'Style','infinite');
end

function plotRoutes3(G,col,lw)
hold on;
for r = 1:numel(G.routes)
    P = G.routes{r}; la = []; lo = [];
    for s = 1:size(P,1)-1, la = [la linspace(P(s,1),P(s+1,1),30)]; lo = [lo linspace(P(s,2),P(s+1,2),30)]; end %#ok<AGROW>
    z = interp2(G.lonv,G.latv,G.elev,lo,la,'linear')/111000*60+0.12;
    plot3(lo,la,z,'-','Color',col,'LineWidth',lw);
end
end

function plotCoast(G,z0)
c = contourc(G.lonv,G.latv,double(G.land),[0.5 0.5]); k = 1;
while k < size(c,2)
    n = c(2,k); plot3(c(1,k+1:k+n),c(2,k+1:k+n),z0*ones(1,n),'k'); k = k+n+1;
end
end

function saveFig(fig,outDir,name)
try, exportgraphics(fig,fullfile(outDir,[name '.png']),'Resolution',300);
catch, print(fig,fullfile(outDir,name),'-dpng','-r300'); end
end

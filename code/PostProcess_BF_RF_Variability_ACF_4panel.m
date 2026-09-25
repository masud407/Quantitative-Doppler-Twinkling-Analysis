% =========================================================
% BF RF VARIABILITY ANALYSIS — ACF + FFT for BOTH
% Pulse-wise (middle 12 of 14) AND Frame-wise (n=20)
%
% Data: fSIG_rf_all → y × x × pulse(14) × frame(20) × rep(5) × sample(5)
%                      1    2      3           4          5        6
%
% Figure layout (4 rows × 2 cols):
%   Row 1: (a) Pulse-wise trajectory       (b) Pulse-wise ACF
%   Row 2: (c) Pulse-wise FFT/PSD          (d) Pulse-wise dominant freq bar
%   Row 3: (e) Frame-wise trajectory       (f) Frame-wise ACF
%   Row 4: (g) Frame-wise FFT/PSD          (h) Frame-wise dominant freq bar
% =========================================================

clear; clc; close all;

% =========================================================
% LOAD DATA
% =========================================================

LPMMA = load('Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\FOC\LPMMA\AverageResults2\LTwinkling_AllSamples_Stacked.mat');
PMMA  = load('Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\FOC\PMMA\AverageResults2\Twinkling_AllSamples_Stacked.mat');
Metal = load('Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\FOC\Metal_Retake\AverageResults2\Metal_Retake_AllSamples_Stacked.mat');

% =========================================================
% EXTRACT BF RF MAGNITUDE
% Dimensions: y(1) x(2) pulse(3) frame(4) rep(5) sample(6)
% =========================================================

PMMA_RF  = abs(PMMA.fSIG_rf_all);
LPMMA_RF = abs(LPMMA.fSIG_rf_all);
Metal_RF = abs(Metal.fSIG_rf_all);

% =========================================================
% SETTINGS
% =========================================================

pulseIdx    = 2:13;           % middle 12 pulses (drop boundary filter edges)
nPulsesUsed = numel(pulseIdx);
nFrames     = 20;
nSamples    = 5;

fontName  = 'Arial';
fontSize  = 12;
labelSize = 14;

cPMMA  = [0.16, 0.47, 0.84];
cLPMMA = [0.07, 0.69, 0.48];
cMetal = [0.89, 0.29, 0.28];
alphaIndiv = 0.18;

outDir = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\FOC\BF_RF_Variability';
if ~exist(outDir,'dir'); mkdir(outDir); end

% =========================================================
% HELPER FUNCTIONS
% =========================================================

cv = @(x) 100 * std(x) / mean(x);

% =========================================================
% PULSE-WISE: mean over y(1) x(2) frame(4) rep(5) → 12 × 5
% =========================================================

PMMA_pulse  = squeeze(mean(mean(mean(mean(PMMA_RF(:,:,pulseIdx,:,:,:),  1),2),4),5));
LPMMA_pulse = squeeze(mean(mean(mean(mean(LPMMA_RF(:,:,pulseIdx,:,:,:), 1),2),4),5));
Metal_pulse = squeeze(mean(mean(mean(mean(Metal_RF(:,:,pulseIdx,:,:,:),  1),2),4),5));

PMMA_pulse_mean  = mean(PMMA_pulse,  2);
LPMMA_pulse_mean = mean(LPMMA_pulse, 2);
Metal_pulse_mean = mean(Metal_pulse, 2);

% =========================================================
% FRAME-WISE: mean over y(1) x(2) pulse(3) rep(5) → 20 × 5
% =========================================================

PMMA_frame  = squeeze(mean(mean(mean(mean(PMMA_RF,  1),2),3),5));
LPMMA_frame = squeeze(mean(mean(mean(mean(LPMMA_RF, 1),2),3),5));
Metal_frame = squeeze(mean(mean(mean(mean(Metal_RF,  1),2),3),5));

PMMA_frame_mean  = mean(PMMA_frame,  2);
LPMMA_frame_mean = mean(LPMMA_frame, 2);
Metal_frame_mean = mean(Metal_frame, 2);

% =========================================================
% ACF — PULSE-WISE (lags 0 … nPulsesUsed-1)
% =========================================================

maxLagP  = nPulsesUsed - 1;
lagAxisP = (0:maxLagP)';

acfP_PMMA  = zeros(maxLagP+1, nSamples);
acfP_LPMMA = zeros(maxLagP+1, nSamples);
acfP_Metal = zeros(maxLagP+1, nSamples);

for s = 1:nSamples
    [r,lg] = xcorr(PMMA_pulse(:,s)  - mean(PMMA_pulse(:,s)),  maxLagP, 'normalized');
    acfP_PMMA(:,s)  = r(lg>=0);
    [r,lg] = xcorr(LPMMA_pulse(:,s) - mean(LPMMA_pulse(:,s)), maxLagP, 'normalized');
    acfP_LPMMA(:,s) = r(lg>=0);
    [r,lg] = xcorr(Metal_pulse(:,s)  - mean(Metal_pulse(:,s)),  maxLagP, 'normalized');
    acfP_Metal(:,s) = r(lg>=0);
end

acfP_PMMA_mean  = mean(acfP_PMMA,  2);
acfP_LPMMA_mean = mean(acfP_LPMMA, 2);
acfP_Metal_mean = mean(acfP_Metal, 2);
acf_confP = 1.96 / sqrt(nPulsesUsed);

% =========================================================
% ACF — FRAME-WISE (lags 0 … nFrames-1)
% =========================================================

maxLagF  = nFrames - 1;
lagAxisF = (0:maxLagF)';

acfF_PMMA  = zeros(maxLagF+1, nSamples);
acfF_LPMMA = zeros(maxLagF+1, nSamples);
acfF_Metal = zeros(maxLagF+1, nSamples);

for s = 1:nSamples
    [r,lg] = xcorr(PMMA_frame(:,s)  - mean(PMMA_frame(:,s)),  maxLagF, 'normalized');
    acfF_PMMA(:,s)  = r(lg>=0);
    [r,lg] = xcorr(LPMMA_frame(:,s) - mean(LPMMA_frame(:,s)), maxLagF, 'normalized');
    acfF_LPMMA(:,s) = r(lg>=0);
    [r,lg] = xcorr(Metal_frame(:,s)  - mean(Metal_frame(:,s)),  maxLagF, 'normalized');
    acfF_Metal(:,s) = r(lg>=0);
end

acfF_PMMA_mean  = mean(acfF_PMMA,  2);
acfF_LPMMA_mean = mean(acfF_LPMMA, 2);
acfF_Metal_mean = mean(acfF_Metal, 2);
acf_confF = 1.96 / sqrt(nFrames);

% =========================================================
% FOUR-PANEL FIGURE
%   (a) Pulse-wise trajectory
%   (b) Pulse-wise ACF
%   (c) Frame-wise trajectory
%   (d) Frame-wise ACF
% =========================================================

fig = figure('Color','w', ...
    'Units','inches','Position',[0.5 0.5 12 9], ...
    'PaperUnits','inches','PaperSize',[12 9], ...
    'PaperPosition',[0 0 12 9]);

panelLabels = {'(a)','(b)','(c)','(d)'};
pulses = (1:nPulsesUsed)';
frames = (1:nFrames)';

% Manual positions [left bottom width height]
lm = 0.09;
rm = 0.56;
w  = 0.38;
h  = 0.35;
bTop = 0.57;
bBot = 0.11;

axPos = [ ...
    lm bTop w h; ...   % (a) pulse trajectory
    rm bTop w h; ...   % (b) pulse ACF
    lm bBot w h; ...   % (c) frame trajectory
    rm bBot w h  ...   % (d) frame ACF
    ];

% =========================================================
% (a) PULSE-WISE TRAJECTORY
% =========================================================
ax(1) = axes('Position',axPos(1,:));
hold(ax(1),'on');
for s = 1:nSamples
    semilogy(ax(1),pulses,PMMA_pulse(:,s), '-','Color',[cPMMA alphaIndiv],'LineWidth',0.8);
    semilogy(ax(1),pulses,LPMMA_pulse(:,s),'-','Color',[cLPMMA alphaIndiv],'LineWidth',0.8);
    semilogy(ax(1),pulses,Metal_pulse(:,s), '-','Color',[cMetal alphaIndiv],'LineWidth',0.8);
end
hP = semilogy(ax(1),pulses,PMMA_pulse_mean, '-','Color',cPMMA,'LineWidth',2.2);
hL = semilogy(ax(1),pulses,LPMMA_pulse_mean,'-','Color',cLPMMA,'LineWidth',2.2);
hM = semilogy(ax(1),pulses,Metal_pulse_mean, '-','Color',cMetal,'LineWidth',2.2);
set(ax(1),'YScale','log','XLim',[1 nPulsesUsed], ...
    'XTick',1:2:nPulsesUsed, ...
    'YLim',[1e1 1e4], ...
    'YTick',[1e1 1e2 1e3 1e4], ...
    'TickDir','out','Box','off', ...
    'FontName',fontName,'FontSize',fontSize,'FontWeight','bold');
grid(ax(1),'on'); ax(1).GridAlpha = 0.15;
xlabel(ax(1),'Pulse number (middle 12)','FontSize',fontSize,'FontWeight','bold','FontName',fontName);
ylabel(ax(1),'Mean BF RF magnitude (a.u.)','FontSize',fontSize,'FontWeight','bold','FontName',fontName);
legend(ax(1),[hP,hL,hM],{'PMMA','LPMMA','Metal'},'Location','best', ...
    'FontSize',10,'FontName',fontName,'Box','on','EdgeColor',[0.7 0.7 0.7]);
text(ax(1),-0.16,1.06,panelLabels{1},'Units','normalized', ...
    'FontName','Times New Roman','FontSize',labelSize,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom');

% =========================================================
% (b) PULSE-WISE ACF
% =========================================================
ax(2) = axes('Position',axPos(2,:));
hold(ax(2),'on');
fill(ax(2),[lagAxisP; flipud(lagAxisP)], ...
    [repmat(acf_confP,maxLagP+1,1); repmat(-acf_confP,maxLagP+1,1)], ...
    [0.85 0.85 0.85],'EdgeColor','none','FaceAlpha',0.50);
for s = 1:nSamples
    plot(ax(2),lagAxisP,acfP_PMMA(:,s), '-','Color',[cPMMA alphaIndiv],'LineWidth',0.8);
    plot(ax(2),lagAxisP,acfP_LPMMA(:,s),'-','Color',[cLPMMA alphaIndiv],'LineWidth',0.8);
    plot(ax(2),lagAxisP,acfP_Metal(:,s), '-','Color',[cMetal alphaIndiv],'LineWidth',0.8);
end
hPa = plot(ax(2),lagAxisP,acfP_PMMA_mean, '-','Color',cPMMA,'LineWidth',2.2);
hLa = plot(ax(2),lagAxisP,acfP_LPMMA_mean,'-','Color',cLPMMA,'LineWidth',2.2);
hMa = plot(ax(2),lagAxisP,acfP_Metal_mean, '-','Color',cMetal,'LineWidth',2.2);
yline(ax(2), acf_confP,'--k','LineWidth',0.9,'Alpha',0.5);
yline(ax(2),-acf_confP,'--k','LineWidth',0.9,'Alpha',0.5);
set(ax(2),'XLim',[0 maxLagP],'XTick',0:2:maxLagP, ...
    'YLim',[-0.55 1], ...
    'TickDir','out','Box','off', ...
    'FontName',fontName,'FontSize',fontSize,'FontWeight','bold');
grid(ax(2),'on'); ax(2).GridAlpha = 0.15;
xlabel(ax(2),'Lag (pulses)','FontSize',fontSize,'FontWeight','bold','FontName',fontName);
ylabel(ax(2),'Normalised ACF','FontSize',fontSize,'FontWeight','bold','FontName',fontName);
legend(ax(2),[hPa,hLa,hMa],{'PMMA','LPMMA','Metal'},'Location','best', ...
    'FontSize',10,'FontName',fontName,'Box','on','EdgeColor',[0.7 0.7 0.7]);
text(ax(2),-0.16,1.06,panelLabels{2},'Units','normalized', ...
    'FontName','Times New Roman','FontSize',labelSize,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom');

% =========================================================
% (c) FRAME-WISE TRAJECTORY
% =========================================================
ax(3) = axes('Position',axPos(3,:));
hold(ax(3),'on');
for s = 1:nSamples
    semilogy(ax(3),frames,PMMA_frame(:,s), '-','Color',[cPMMA alphaIndiv],'LineWidth',0.8);
    semilogy(ax(3),frames,LPMMA_frame(:,s),'-','Color',[cLPMMA alphaIndiv],'LineWidth',0.8);
    semilogy(ax(3),frames,Metal_frame(:,s), '-','Color',[cMetal alphaIndiv],'LineWidth',0.8);
end
hP3 = semilogy(ax(3),frames,PMMA_frame_mean, '-','Color',cPMMA,'LineWidth',2.2);
hL3 = semilogy(ax(3),frames,LPMMA_frame_mean,'-','Color',cLPMMA,'LineWidth',2.2);
hM3 = semilogy(ax(3),frames,Metal_frame_mean, '-','Color',cMetal,'LineWidth',2.2);
set(ax(3),'YScale','log','XLim',[1 nFrames], ...
    'XTick',[1 5 10 15 20], ...
    'YLim',[1e1 1e4], ...
    'YTick',[1e1 1e2 1e3 1e4], ...
    'TickDir','out','Box','off', ...
    'FontName',fontName,'FontSize',fontSize,'FontWeight','bold');
grid(ax(3),'on'); ax(3).GridAlpha = 0.15;
xlabel(ax(3),'Frame number','FontSize',fontSize,'FontWeight','bold','FontName',fontName);
ylabel(ax(3),'Mean BF RF magnitude (a.u.)','FontSize',fontSize,'FontWeight','bold','FontName',fontName);
legend(ax(3),[hP3,hL3,hM3],{'PMMA','LPMMA','Metal'},'Location','best', ...
    'FontSize',10,'FontName',fontName,'Box','on','EdgeColor',[0.7 0.7 0.7]);
text(ax(3),-0.16,1.06,panelLabels{3},'Units','normalized', ...
    'FontName','Times New Roman','FontSize',labelSize,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom');

% =========================================================
% (d) FRAME-WISE ACF
% =========================================================
ax(4) = axes('Position',axPos(4,:));
hold(ax(4),'on');
fill(ax(4),[lagAxisF; flipud(lagAxisF)], ...
    [repmat(acf_confF,maxLagF+1,1); repmat(-acf_confF,maxLagF+1,1)], ...
    [0.85 0.85 0.85],'EdgeColor','none','FaceAlpha',0.50);
for s = 1:nSamples
    plot(ax(4),lagAxisF,acfF_PMMA(:,s), '-','Color',[cPMMA alphaIndiv],'LineWidth',0.8);
    plot(ax(4),lagAxisF,acfF_LPMMA(:,s),'-','Color',[cLPMMA alphaIndiv],'LineWidth',0.8);
    plot(ax(4),lagAxisF,acfF_Metal(:,s), '-','Color',[cMetal alphaIndiv],'LineWidth',0.8);
end
hPf = plot(ax(4),lagAxisF,acfF_PMMA_mean, '-','Color',cPMMA,'LineWidth',2.2);
hLf = plot(ax(4),lagAxisF,acfF_LPMMA_mean,'-','Color',cLPMMA,'LineWidth',2.2);
hMf = plot(ax(4),lagAxisF,acfF_Metal_mean, '-','Color',cMetal,'LineWidth',2.2);
yline(ax(4), acf_confF,'--k','LineWidth',0.9,'Alpha',0.5);
yline(ax(4),-acf_confF,'--k','LineWidth',0.9,'Alpha',0.5);
set(ax(4),'XLim',[0 maxLagF],'XTick',0:5:maxLagF, ...
    'YLim',[-0.55 1], ...
    'TickDir','out','Box','off', ...
    'FontName',fontName,'FontSize',fontSize,'FontWeight','bold');
grid(ax(4),'on'); ax(4).GridAlpha = 0.15;
xlabel(ax(4),'Lag (frames)','FontSize',fontSize,'FontWeight','bold','FontName',fontName);
ylabel(ax(4),'Normalised ACF','FontSize',fontSize,'FontWeight','bold','FontName',fontName);
legend(ax(4),[hPf,hLf,hMf],{'PMMA','LPMMA','Metal'},'Location','best', ...
    'FontSize',10,'FontName',fontName,'Box','on','EdgeColor',[0.7 0.7 0.7]);
text(ax(4),-0.16,1.06,panelLabels{4},'Units','normalized', ...
    'FontName','Times New Roman','FontSize',labelSize,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom');

% =========================================================
% SAVE
% =========================================================

exportgraphics(fig, fullfile(outDir,'BF_RF_Pulse_Frame_ACF_4panel.png'), ...
    'Resolution',600,'BackgroundColor','white');
exportgraphics(fig, fullfile(outDir,'BF_RF_Pulse_Frame_ACF_4panel.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
savefig(fig, fullfile(outDir,'BF_RF_Pulse_Frame_ACF_4panel.fig'));

fprintf('\n4-panel figure saved to: %s\n', outDir);

clear; clc; close all;

% =========================================================================
%  USER PATHS
% =========================================================================
pmma_dir  = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\PW\PMMA\Each1\X1-X2_Investigation';
lpmma_dir = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\PW\LPMMA\Each1\X1-X2_Investigation';
metal_dir = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\PW\Metal_Retake\Each1\X1-X2_Investigation';
% 
% pmma_dir  = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\PW\PMMA\All\X1-X2_Investigation';
% lpmma_dir = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\PW\LPMMA\All\X1-X2_Investigation';
% metal_dir = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\PW\Metal_Retake\All\X1-X2_Investigation';

% metal_dir  = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L7-4\PW\PMMA\S4\6MHz\S0\X1-X2_Investigation';
% lpmma_dir = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L7-4\PW\LPMMA\S4\6MHz\S0\X1-X2_Investigation';
%  pmma_dir = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L7-4\PW\Metal_Retake\S5\6 MHz\S0\X1-X2_Investigation';
%  % metal_dir = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L7-4\PW\Metal_Retake\Each1\X1-X2_Investigation';
 out_root = 'Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\PW\Publication_Figures';
if ~exist(out_root,'dir'), mkdir(out_root); end

% =========================================================================
%  APPEARANCE
%  KEY FIXES:
%  - N_PTS reduced to 800  → fast rendering
%  - ALPHA_SC raised to 0.55 → visible colors
%  - MS_SC = 3              → smaller dots, less overlap
% =========================================================================
% FIG_2D_W  = 20;   % cm double-column
% FIG_2D_H  =  8;
% FIG_3D_W  = 15.0;
% FIG_3D_H  =  10.0;
FIG_2D_W  = 40;   % cm double-column
FIG_2D_H  =  12;
FIG_3D_W  = 28.0;
FIG_3D_H  =  19.0;
DPI       = 600;

% FS_PANEL  = 10;
% FS_LABEL  =  9;
% FS_TICK   =  8;
% FS_LEGEND =  8;
FS_PANEL  = 18;
FS_LABEL  =  18;
FS_TICK   =  18;
FS_LEGEND =  16;
LW_AXIS   =  0.8;
MS_SC     =  5;      % ← smaller markers
MS_MEAN   = 9;
ALPHA_SC  =  0.55;   % ← raised: colors now clearly visible
N_PTS     =  800;    % ← reduced: fast render, no slowdown
N_PTS_3D  =  600;

MAT_LABELS = {'PMMA','LPMMA','Metal'};
MAT_DIRS   = {pmma_dir, lpmma_dir, metal_dir};

% Colorblind-safe: blue / orange / teal
% C = [0.09  0.47  0.71;
%      0.93  0.54  0.00;
%      0.17  0.63  0.60];
C = [0.00 0.45 0.74;   % PMMA  — blue
          0.20 0.65 0.20;   % LPMMA — green
          0.80 0.15 0.15];  % Metal  — red

C_EDGE = [0 0 0;   % PMMA — black
          0 0 0;   % LPMMA — black
          0 0 0];  % Metal — black

% ── at the top with color definitions ──────────────────────
% C_MEAN_FILL = [0.85 0.10 0.10];   % bright red — all 3 means
% C_MEAN_EDGE = [0.50 0.00 0.00];   % dark red edge
% C_MEAN_HALO = [1.00 1.00 1.00];   % white halo

MSTYLE = {'o','s','^'};

N3D_BINS  = 25;
pd3_edges = linspace(3,   11,   N3D_BINS+1);
x2_3edges = linspace(0,   1.60, N3D_BINS+1);
x1_3edges = linspace(0,   0.65, N3D_BINS+1);
pairs = {[1 2],[1 3],[2 3]};

% =========================================================================
%  STEP 1 — Load data
% =========================================================================
fprintf('Loading data ...\n');
file_lists = cell(3,1);
for m = 1:3
    fl = dir(fullfile(MAT_DIRS{m},'*.mat'));
    if isempty(fl), error('No .mat files in: %s', MAT_DIRS{m}); end
    [~,si] = sort({fl.name});
    file_lists{m} = fl(si);
end
n_samples = min(cellfun(@length, file_lists));
fprintf('  %d samples per material.\n\n', n_samples);

DATA = cell(3, n_samples);
POOL = cell(3,1);
for m = 1:3
    POOL{m} = struct('x1',[],'x2',[],'logPD',[],'x1m',[],'x2m',[],'PD',[]);
end

for m = 1:3
    for s = 1:n_samples
        raw  = load(fullfile(MAT_DIRS{m}, file_lists{m}(s).name));
        PD   = raw.denom_all(:);
        x1   = abs(raw.num_im_all(:) ./ raw.denom_all(:));
        x2   = 1 - raw.num_re_all(:) ./ raw.denom_all(:);
        mask = PD > 0;
        DATA{m,s} = struct('x1',x1,'x2',x2,'PD',PD, ...
                           'logPD',log10(PD(mask)), ...
                           'x1m',x1(mask),'x2m',x2(mask));
        POOL{m}.x1    = [POOL{m}.x1;    x1];
        POOL{m}.x2    = [POOL{m}.x2;    x2];
        POOL{m}.PD    = [POOL{m}.PD;    PD];
        POOL{m}.logPD = [POOL{m}.logPD; log10(PD(mask))];
        POOL{m}.x1m   = [POOL{m}.x1m;  x1(mask)];
        POOL{m}.x2m   = [POOL{m}.x2m;  x2(mask)];
        fprintf('  Loaded [%s] sample %d\n', MAT_LABELS{m}, s);
    end
end
fprintf('\nAll data loaded.\n\n');

% =========================================================================
%  STEP 2 — BC3D table
% =========================================================================
fprintf('Computing BC3D ...\n');
row_labels={}; x1_str={}; x2_str={}; pd_str={};
bc_12=[]; bc_13=[]; bc_23=[];

for s = 1:n_samples
    Ds = DATA(:,s);
    bc3s = zeros(1,3);
    for p = 1:3
        i=pairs{p}(1); j=pairs{p}(2);
        mi=Ds{i}.PD>0; mj=Ds{j}.PD>0;
        P3=hist3d(log10(Ds{i}.PD(mi)),Ds{i}.x2(mi),Ds{i}.x1(mi),pd3_edges,x2_3edges,x1_3edges);
        Q3=hist3d(log10(Ds{j}.PD(mj)),Ds{j}.x2(mj),Ds{j}.x1(mj),pd3_edges,x2_3edges,x1_3edges);
        P3=P3/sum(P3(:)); Q3=Q3/sum(Q3(:));
        bc3s(p)=sum(sqrt(P3(:).*Q3(:)));
    end
    for m = 1:3
        d=Ds{m};
        row_labels{end+1} = sprintf('%s — Sample %d',MAT_LABELS{m},s); %#ok
        x1_str{end+1}     = sprintf('%.3f ± %.3f',mean(d.x1),   std(d.x1));    %#ok
        x2_str{end+1}     = sprintf('%.3f ± %.3f',mean(d.x2),   std(d.x2));    %#ok
        pd_str{end+1}     = sprintf('%.3f ± %.3f',mean(d.logPD),std(d.logPD)); %#ok
        bc_12(end+1)=bc3s(1); bc_13(end+1)=bc3s(2); bc_23(end+1)=bc3s(3);     %#ok
    end
    fprintf('  Sample %d BC3D: [%.4f | %.4f | %.4f]\n',s,bc3s(1),bc3s(2),bc3s(3));
end
for m = 1:3
    mat_rows = m:3:(n_samples*3);
    d=POOL{m};
    row_labels{end+1} = sprintf('**%s — Average**',MAT_LABELS{m}); %#ok
    x1_str{end+1}     = sprintf('%.3f ± %.3f',mean(d.x1),   std(d.x1));    %#ok
    x2_str{end+1}     = sprintf('%.3f ± %.3f',mean(d.x2),   std(d.x2));    %#ok
    pd_str{end+1}     = sprintf('%.3f ± %.3f',mean(d.logPD),std(d.logPD)); %#ok
    bc_12(end+1)=mean(bc_12(mat_rows)); %#ok
    bc_13(end+1)=mean(bc_13(mat_rows)); %#ok
    bc_23(end+1)=mean(bc_23(mat_rows)); %#ok
end

T = table(row_labels(:),x1_str(:),x2_str(:),pd_str(:), ...
          bc_12(:),bc_13(:),bc_23(:), ...
    'VariableNames',{'Sample','x1_mean_SD','x2_mean_SD','log10PD_mean_SD', ...
                     'BC3D_PMMA_LPMMA','BC3D_PMMA_Metal','BC3D_LPMMA_Metal'});
writetable(T, fullfile(out_root,'Table_PerSample_Summary.csv'));
writetable(T, fullfile(out_root,'Table_PerSample_Summary.xlsx'));
fprintf('Saved table.\n\n');

% =========================================================================
%  FIGURE 1 — 1×3 panels (a)(b)(c)
%  Render strategy: plot() instead of scatter() for speed,
%  then overlay mean markers with scatter()
% =========================================================================
fprintf('Generating 2-D figure (a)(b)(c) ...\n');

rng(0);

fig2d = figure('Color','w', ...
               'Units','centimeters', ...
               'Position',[2 2 FIG_2D_W FIG_2D_H], ...
               'PaperUnits','centimeters', ...
               'PaperSize',[FIG_2D_W FIG_2D_H], ...
               'PaperPosition',[0 0 FIG_2D_W FIG_2D_H]);

% Panel specs
%  col 1-5: xlabel ylabel xlim ylim label
%  col 6-7: xdata-getter  ydata-getter  (applied to POOL{m})
pXL   = {'x_2  (spectral broadening)',     'log_{10}(PD)',  'log_{10}(PD)'};
pYL   = {'x_1  (mean Doppler frequency)',      'x_1  (mean Doppler frequency)', 'x_2  (spectral broadening)'};
pXLIM = {[0 1.2],  [3 11],   [3 11]};
pYLIM = {[0 0.30], [0 0.65], [0 1.60]};
pLAB  = {'(a)','(b)','(c)'};

getX  = { @(d) d.x2,    @(d) d.logPD, @(d) d.logPD };
getY  = { @(d) d.x1,    @(d) d.x1m,   @(d) d.x2m   };

lm=0.085; rm=0.015; bm=0.19; tm=0.06; gap=0.065;
pw = (1 - lm - rm - 2*gap) / 3;
ph = 1 - bm - tm;
col_l = lm + (0:2).*(pw+gap);

for p = 1:3
    ax = axes('Parent',fig2d, ...                    %#ok<LAXES>
              'Position',[col_l(p) bm pw ph]);
    hold(ax,'on');

    h_leg = gobjects(1,3);

    for m = 1:3
        xall = getX{p}(POOL{m});
        yall = getY{p}(POOL{m});

        % subsample
        np  = min(N_PTS, numel(xall));
        idx = randperm(numel(xall), np);
        xs  = xall(idx);
        ys  = yall(idx);

        % ── USE plot() NOT scatter() → 10-20× faster ─────────────────
        h_leg(m) = plot(ax, xs, ys, '.', ...
                        'Color',     [C(m,:) ALPHA_SC], ...
                        'MarkerSize', MS_SC+2, ...
                        'DisplayName', MAT_LABELS{m});
    end

    % ── Mean markers — white halo then filled symbol ──────────────────
    for m = 1:3
        mx = mean(getX{p}(POOL{m}));
        my = mean(getY{p}(POOL{m}));

        % outer white ring for contrast
        plot(ax, mx, my, MSTYLE{m}, ...
             'Color',           'w',   ...
             'MarkerFaceColor', 'w',   ...
             'MarkerSize',      MS_MEAN+4, ...
             'LineWidth',       1,     ...
             'HandleVisibility','off');

        % filled colored marker with dark edge
        plot(ax, mx, my, MSTYLE{m}, ...
             'Color',           C_EDGE(m,:), ...
             'MarkerFaceColor', C(m,:),       ...
             'MarkerSize',      MS_MEAN,      ...
             'LineWidth',       1,          ...
             'HandleVisibility','off');
    end

    % ── axes style ────────────────────────────────────────────────────
    xlim(ax, pXLIM{p});
    ylim(ax, pYLIM{p});

    xlabel(ax, pXL{p}, 'FontSize',FS_LABEL, 'FontWeight','bold', ...
               'FontName','Arial', 'Interpreter','tex');
    ylabel(ax, pYL{p}, 'FontSize',FS_LABEL, 'FontWeight','bold', ...
               'FontName','Arial', 'Interpreter','tex');

    set(ax, 'FontSize',      FS_TICK,  ...
            'FontName',      'Arial',  ...
            'LineWidth',     LW_AXIS,  ...
            'FontWeight',     'bold',         ...
            'Box',           'on',     ...
            'TickDir',       'out',    ...
            'XGrid',         'on',     ...
            'YGrid',         'on',     ...
            'GridAlpha',     0.20,     ...
            'GridColor',     [0.5 0.5 0.5], ...
            'GridLineStyle', ':');

    % panel label
    xl=pXLIM{p}; yl=pYLIM{p};
    % text(ax, xl(1)+diff(xl)*0.03, yl(2)-diff(yl)*0.04, pLAB{p}, ...
    %      'FontSize',  FS_PANEL+1, 'FontWeight','bold', ...
    %      'FontName',  'Arial',    'Interpreter','none', ...
    %      'HorizontalAlignment','left','VerticalAlignment','top');
    % text(ax, -0.10, 1.03, pLAB{p}, ...
    % 'Units','normalized', ...
    % 'FontSize',FS_PANEL+1, ...
    % 'FontWeight','bold', ...
    % 'FontName','Arial', ...
    % 'HorizontalAlignment','left', ...
    % 'VerticalAlignment','bottom', ...
    % 'Clipping','off');

    % legend on panel (a) only
    if p == 1
        lg = legend(ax, h_leg, MAT_LABELS, ...
                    'FontSize',  FS_LEGEND, ...
                    'FontName',  'Arial',   ...
                    'Location',  'northeast', ...
                    'Box',       'on',      ...
                    'EdgeColor', [0.65 0.65 0.65]);
        lg.ItemTokenSize = [8 8];
        % text(ax, xl(2)-diff(xl)*0.03, yl(1)+diff(yl)*0.03, ...
        %      'Large symbols = group mean', ...
        %      'FontSize', FS_TICK-0.5, 'FontName','Arial', ...
        %      'Color',[0.4 0.4 0.4], ...
        %      'HorizontalAlignment','right','VerticalAlignment','bottom');
    end

    fprintf('  Panel %s done.\n', pLAB{p});
end

fprintf('Saving 2-D figure ...\n');
exportgraphics(fig2d, fullfile(out_root,'Fig_Scatter_2D_abc.png'), ...
    'Resolution',DPI,'BackgroundColor','white');
exportgraphics(fig2d, fullfile(out_root,'Fig_Scatter_2D_abc.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
savefig(fig2d, fullfile(out_root,'Fig_Scatter_2D_abc.fig'));
fprintf('  Saved: Fig_Scatter_2D_abc\n\n');

% =========================================================================
%  FIGURE 2 — Standalone 3-D scatter (d)
% =========================================================================
fprintf('Generating 3-D figure (d) ...\n');

fig3d = figure('Color','w', ...
               'Units','centimeters', ...
               'Position',[14 2 FIG_3D_W FIG_3D_H], ...
               'PaperUnits','centimeters', ...
               'PaperSize',[FIG_3D_W FIG_3D_H], ...
               'PaperPosition',[0 0 FIG_3D_W FIG_3D_H]);

% ax3d = axes('Parent',fig3d,'Position',[0.10 0.08 0.85 0.84]);
ax3d = axes('Parent',fig3d, ...
            'Position',[0.1 0.1 0.76 0.76]);
hold(ax3d,'on');

h3 = gobjects(1,3);
for m = 1:3
    np  = min(N_PTS_3D, numel(POOL{m}.logPD));
    idx = randperm(numel(POOL{m}.logPD), np);

    % plot3() instead of scatter3() — much faster
    h3(m) = plot3(ax3d, ...
                  POOL{m}.logPD(idx), ...
                  POOL{m}.x2m(idx),   ...
                  POOL{m}.x1m(idx),   ...
                  '.', ...
                  'Color',      [C(m,:) ALPHA_SC], ...
                  'MarkerSize', MS_SC+3,            ...
                  'DisplayName', MAT_LABELS{m});
end

% % mean markers with white halo
% for m = 1:3
%     mx = mean(POOL{m}.logPD);
%     my = mean(POOL{m}.x2m);
%     mz = mean(POOL{m}.x1m);
% 
%     % white halo
%     % plot3(ax3d, mx, my, mz, MSTYLE{m}, ...
%     %       'Color','w','MarkerFaceColor','w', ...
%     %       'MarkerSize',MS_MEAN+7,'LineWidth',3, ...
%     %       'HandleVisibility','off');
% 
%     plot3(ax3d, mx, my, mz, MSTYLE{m}, ...
%     'Color',            [1 1 1],   ...   % ← white edge
%     'MarkerFaceColor',  C(m,:),    ...   % ← same color as scatter
%     'MarkerEdgeColor',  [0 0 0],   ...   % ← white edge
%     'MarkerSize',       MS_MEAN+4, ...
%     'LineWidth',        3,       ...
%     'HandleVisibility', 'off');
% 
% 
% % filled marker — change C(m,:) to C_MEAN_FILL and C_EDGE(m,:) to C_MEAN_EDGE
% % plot3(ax3d, mx, my, mz, MSTYLE{m}, ...
% %     'Color',           [0.50 0.00 0.00], ...   % ← dark red edge
% %     'MarkerFaceColor', [0.85 0.10 0.10], ...   % ← bright red fill
% %     'MarkerSize',      MS_MEAN+2,        ...
% %     'LineWidth',       1.5,              ...
% %     'HandleVisibility','off');
% 
% % text label
% % text(ax3d, mx, my, mz+0.025, MAT_LABELS{m}, ...
% %     'FontSize',  FS_LEGEND+2,      ...
% %     'FontWeight','bold',           ...
% %     'FontName',  'Arial',          ...
% %     'Color',     [0.50 0.00 0.00], ...   % ← dark red label
% %     'HorizontalAlignment','center');
% 
% % text(ax3d, mx, my, mz+0.025, MAT_LABELS{m}, ...
% %     'FontSize',  FS_LEGEND+2, ...
% %     'FontWeight','bold',      ...
% %     'FontName',  'Arial',     ...
% %     'Color',     C(m,:),      ...   % ← matches scatter color
% %     'HorizontalAlignment','center');
% end
for m = 1:3

    mx = mean(POOL{m}.logPD);
    my = mean(POOL{m}.x2m);
    mz = mean(POOL{m}.x1m);

    % ---------------------------------------------------------
    % White halo behind mean marker
    % ---------------------------------------------------------
    plot3(ax3d, mx, my, mz, MSTYLE{m}, ...
        'Color',            [1 1 1], ...
        'MarkerFaceColor',  [1 1 1], ...
        'MarkerEdgeColor',  [1 1 1], ...
        'MarkerSize',       MS_MEAN + 9, ...
        'LineWidth',        4, ...
        'HandleVisibility', 'off');

    % ---------------------------------------------------------
    % Mean marker
    % ---------------------------------------------------------
    plot3(ax3d, mx, my, mz, MSTYLE{m}, ...
        'Color',            [0 0 0], ...       % BLACK EDGE
        'MarkerFaceColor',  C(m,:), ...         % material color
        'MarkerEdgeColor',  [0 0 0], ...       % BLACK EDGE
        'MarkerSize',       MS_MEAN + 5, ...
        'LineWidth',        2.5, ...
        'HandleVisibility', 'off');
end
% xlabel(ax3d,'log_{10}(PD)',               'FontSize',FS_LABEL+1,'FontWeight','bold','FontName','Arial','Interpreter','tex');
% ylabel(ax3d,'x_2  (spectral broadening)', 'FontSize',FS_LABEL+1,'FontWeight','bold','FontName','Arial','Interpreter','tex');
% ylab = ylabel(ax3d,'x_2  (spectral broadening)', ...
%     'FontSize',FS_LABEL, ...
%     'FontWeight','bold', ...
%     'FontName','Arial', ...
%     'Interpreter','tex');
% 
% ylab.Rotation = 15;

% X-axis label
xlabel(ax3d,'log_{10}(PD)', ...
    'FontSize',FS_LABEL+1, ...
    'FontWeight','bold', ...
    'FontName','Arial', ...
    'Interpreter','tex');

% Y-axis label: manually position and rotate
yl = ylabel(ax3d,'x_2  (spectral broadening)', ...
    'FontSize',FS_LABEL+1, ...
    'FontWeight','bold', ...
    'FontName','Arial', ...
    'Interpreter','tex');

% Move it closer to the x2 axis
yl.Position = [11.9, 0.77, 0];

% Rotate to approximately follow the projected x2 axis
yl.Rotation = 22;
zlabel(ax3d,'x_1  (mean Doppler frequency)',  'FontSize',FS_LABEL+1,'FontWeight','bold','FontName','Arial','Interpreter','tex');

set(ax3d,'FontSize',  FS_TICK+1, ...
         'FontName',  'Arial',   ...
         'FontWeight',     'bold',         ...
         'LineWidth', LW_AXIS,   ...
         'Box',       'on',      ...
         'XGrid','on','YGrid','on','ZGrid','on', ...
         'GridAlpha', 0.20,      ...
         'GridColor', [0.5 0.5 0.5], ...
         'GridLineStyle',':');

view(ax3d, 38, 26);

lg3d = legend(ax3d, h3, MAT_LABELS, ...
              'FontSize',  FS_LEGEND+1, ...
              'FontName',  'Arial',     ...
              'Location',  'best', ...
              'Box',       'on',        ...
              'EdgeColor', [0.65 0.65 0.65]);
lg3d.ItemTokenSize = [10 10];

% title(ax3d,'(d)','FontSize',FS_PANEL+2,'FontWeight','bold', ...
%      'FontName','Arial','HorizontalAlignment','left');

fprintf('Saving 3-D figure ...\n');
exportgraphics(fig3d, fullfile(out_root,'Fig_Scatter_3D_d.png'), ...
    'Resolution',DPI,'BackgroundColor','white');
exportgraphics(fig3d, fullfile(out_root,'Fig_Scatter_3D_d.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
savefig(fig3d, fullfile(out_root,'Fig_Scatter_3D_d.fig'));
fprintf('  Saved: Fig_Scatter_3D_d\n\n');

fprintf('=========================================================\n');
fprintf('  DONE. Outputs in:\n  %s\n', out_root);
fprintf('=========================================================\n');

% =========================================================================
%  LOCAL FUNCTION
% =========================================================================
function H = hist3d(A,B,C,edgesA,edgesB,edgesC)
    nA=numel(edgesA)-1; nB=numel(edgesB)-1; nC=numel(edgesC)-1;
    ia=discretize(A,edgesA); ib=discretize(B,edgesB); ic=discretize(C,edgesC);
    v=~isnan(ia)&~isnan(ib)&~isnan(ic);
    lin=sub2ind([nA nB nC],ia(v),ib(v),ic(v));
    H=reshape(accumarray(lin,1,[nA*nB*nC 1]),[nA nB nC]);
end
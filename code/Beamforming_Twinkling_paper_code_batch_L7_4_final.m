clear all;
close all;
clc;

tic;

addpath('MUST/');

%% =====================================================
% MAIN ROOT FOLDER
%% =====================================================

 % rootFolder ='Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L7-4\Foc\LPMMA\6MHz';
rootFolder ='Y:\Masud\Twinkling_signal_process\Twinkle_paper_Data\L7-4\PW\Metal_Retake\S5\6 MHz';
% 'Y:\URBAN_RENAL_R01\Masud\Twinkling_signal_process\Twinkle_paper_Data\L11-4\PW\LPMMA';


%% =====================================================
% FIXED ROI COORDINATES
%% =====================================================

x_start_user = 50;
x_end_user   = 80;

% y_start_user = 208;
% y_end_user   = 225;

% for case 2
% y_start_user = 208;
% y_end_user   = 232;

% for case 2
% y_start_user = 200;
% y_end_user   = 232;
%for case 1
% y_start_user = 175;  
% y_end_user   = 195;

%for case 2 (more depth ROI)
 % y_start_user = 175;  
 % y_end_user   = 210;

 %for case 2 (more depth ROI)
 y_start_user = 165;  
 y_end_user   = 203;

%% =====================================================
% FIND SAMPLE FOLDERS
%% =====================================================

sampleFolders = dir(fullfile(rootFolder,'S*'));

sampleFolders = sampleFolders([sampleFolders.isdir]);

fprintf('\nFound %d sample folders\n',length(sampleFolders));

%% =====================================================
% LOOP OVER SAMPLE FOLDERS
%% =====================================================

for sIdx = 1:length(sampleFolders)

    sampleName = sampleFolders(sIdx).name;

    samplePath = fullfile(rootFolder,sampleName);

    fprintf('\n=====================================\n');
    fprintf('Processing Sample: %s\n',sampleName);
    fprintf('=====================================\n');

    %% -------------------------------------------------
    % OUTPUT FOLDER
    %% -------------------------------------------------

    outputFolder = ...
        fullfile(samplePath,'ProcessedResults2');  % 2 for bigger ROI_100 HZ WF, 3 is for bigger ROI, EF_200Hz

    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    %% -------------------------------------------------
    % FIND CHANNEL FILES
    %% -------------------------------------------------

    chFiles = dir(fullfile(samplePath,...
        'CHData_B_Doppler_*.mat'));

    fprintf('Found %d channel files\n',length(chFiles));

    %% =================================================
    % LOOP OVER REPETITIONS
    %% =================================================

    for repIdx = 1:length(chFiles)

        close all;

        fprintf('\n-------------------------------------\n');
        fprintf('Processing repetition %d / %d\n', ...
            repIdx,length(chFiles));
        fprintf('-------------------------------------\n');

        %% =============================================
        % LOAD CHANNEL DATA
        %% =============================================

        fileName = chFiles(repIdx).name;

        filePath = fullfile(samplePath,fileName);

        fprintf('Loading:\n%s\n',fileName);

        load(filePath);

        %% -------------------------------------------------
        % BASE NAME
        %% -------------------------------------------------

        [~,baseName,~] = fileparts(fileName);

        resultBase = sprintf('%s_Rep_%02d',...
            sampleName,repIdx);

        %% =============================================
        % PARAMETERS
        %% =============================================

        angles = 7;
        nD = 14;
        f0 = 4e6;  %%
        fs = f0*4;
        NFrames = 20;
         rfB = 2048;
         rfD = 2048;

        %      rfB = 1664;
        % rfD = 1664;

        %% =============================================
        % PARSE RF DATA
        %% =============================================

        rangeB = zeros(angles,2);

        rf_dataB = zeros(rfB,128,angles,NFrames);

        for n = 1:NFrames

            for nn = 1:angles

                rangeB(nn,:) = ...
                    [1+((nn-1)*rfB),nn*rfB];

                rf_dataB(:,:,nn,n) = ...
                    double( ...
                    ChData_B_Doppler( ...
                    rangeB(nn,1):rangeB(nn,2),:,n));

            end
        end

        rangeD = zeros(nD,2);

        rf_dataD = zeros(rfD,128,angles,NFrames);

        for n = 1:NFrames

            for nn = 1:nD

                rangeD(nn,:) = ...
                    [(rfB*7+1)+((nn-1)*rfD), ...
                    (rfB*7)+nn*rfD];

                rf_dataD(:,:,nn,n) = ...
                    double( ...
                    ChData_B_Doppler( ...
                    rangeD(nn,1):rangeD(nn,2),:,n));

            end
        end

        %% =============================================
        % BEAMFORMING PARAMETERS
        %% =============================================

        % fs = 20e6;

        param.fs = fs;

        param.pitch = 0.308e-3;

        param.fc = f0;

        param.c = 1540;

        param.fnumber = 1.5;

        param.TXdelay = zeros(1,128);

        param.radius = Inf;

        dx = param.pitch;

        dz = 0.5*param.pitch;

        x = 0:dx:dx*123;

        x = x - max(x)/2;

        z = 0:dz:dz*512;

        [X,Z] = meshgrid(x,z);

        CF = 2e-07;

        Fs_original = fs;

        Fc = f0;

        param.PRF = 2.5e3;

        bfSIG_rf = zeros(length(z),length(x),14,20);

        bfSIG_IQ = [];

        fSIG = [];

        [bf,af] = butter(4,100/param.PRF/2,'high');

        %% =============================================
        % MAIN PROCESSING
        %% =============================================

        for n = 1:20

            fprintf('Frame %d / 20\n',n);

            SIG_rfD = double(rf_dataD(:,:,:,n));

            param.TXdelay = CF*(TX(8).Delay);

            bfSIG_rf(:,:,:,n) = ...
                das(SIG_rfD,X,Z,...
                param.TXdelay,param);

            bfSIG_IQ_20 = ...
                rf2iq( ...
                bfSIG_rf(:,:,:,n), ...
                Fs_original, ...
                Fc);

            bfSIG_IQ(:,:,:,n) = ...
                bfSIG_IQ_20;

            for n1 = 1:size(bfSIG_IQ,1)

                for n2 = 1:size(bfSIG_IQ,2)

                    temp = ...
                        squeeze( ...
                        bfSIG_IQ(n1,n2,:,n));
                                temp2 = squeeze( bfSIG_rf(n1,n2,:,n));

                    fSIG(n1,n2,:,n) = ...
                      filtfilt(bf,af,temp);
                    fSIG_rf(n1,n2,:,n) =filtfilt(bf,af,temp2);

                end
            end
        end

        %% =============================================
        % DOPPLER PROCESSING
        %% =============================================

        % Fs = 20e6;

        param.fs = fs;

        z = z(1:2:end);

        [X,Z] = meshgrid(x,z);

        for n = 1:20

            Doppler(:,:,n) = ...
                iq2doppler( ...
                bfSIG_IQ(:,:,:,n), ...
                param,[7 7]);

            DopplerWF(:,:,n) = ...
                iq2doppler( ...
                fSIG(:,:,:,n), ...
                param,[7 7]);

            PD(:,:,n) = ...
                sum(abs(bfSIG_IQ(:,:,:,n)).^2,3);

            PDWF(:,:,n) = ...
                sum(abs(fSIG(:,:,:,n)).^2,3);

        end

        %% =============================================
        % SAFE ROI LIMITS
        %% =============================================

        [nz,nx,~] = size(PD);

        x_start = max(1,x_start_user);
        x_end   = min(nx,x_end_user);

        y_start = max(1,y_start_user);
        y_end   = min(nz,y_end_user);

        fprintf('\nROI USED:\n');
        fprintf('x : %d -> %d\n',x_start,x_end);
        fprintf('y : %d -> %d\n',y_start,y_end);

        %% =============================================
        % DOPPLER DISPLAY VIDEO
        %% =============================================

        videoName = fullfile(outputFolder,...
            [resultBase '_DopplerDisplay.avi']);

        v = VideoWriter(videoName);

        open(v);

        figD = figure(99);

        for n = 1:20

            clf;

            set(gcf,'Color','w');

            subplot(221)

            imagesc(Doppler(:,:,n))

            axis image ij tight

            title('Color Doppler')

            colormap dopplermap

            colorbar

            hold on

            rectangle( ...
                'Position',...
                [x_start ...
                 y_start ...
                 x_end-x_start ...
                 y_end-y_start],...
                 'EdgeColor','r',...
                 'LineWidth',2);

            subplot(222)

            imagesc(DopplerWF(:,:,n))

            axis image ij tight

            title('Color Doppler WF')

            colormap dopplermap

            colorbar

            hold on

            rectangle( ...
                'Position',...
                [x_start ...
                 y_start ...
                 x_end-x_start ...
                 y_end-y_start],...
                 'EdgeColor','r',...
                 'LineWidth',2);

            subplot(223)

            imagesc(PD(:,:,n))

            axis image ij tight

            title('Power Doppler')

            colormap hot

            colorbar

            hold on

            rectangle( ...
                'Position',...
                [x_start ...
                 y_start ...
                 x_end-x_start ...
                 y_end-y_start],...
                 'EdgeColor','g',...
                 'LineWidth',2);

            subplot(224)

            imagesc(PDWF(:,:,n))

            axis image ij tight

            title('Power Doppler WF')

            colormap hot

            colorbar

            hold on

            rectangle( ...
                'Position',...
                [x_start ...
                 y_start ...
                 x_end-x_start ...
                 y_end-y_start],...
                 'EdgeColor','g',...
                 'LineWidth',2);

            drawnow;

            frame = getframe(gcf);

            writeVideo(v,frame);

        end

        close(v);

        saveas(figD,...
            fullfile(outputFolder,...
            [resultBase '_DopplerDisplay.png']));

        savefig(figD,...
            fullfile(outputFolder,...
            [resultBase '_DopplerDisplay.fig']));

        %% =============================================
        % ROI EXTRACTION
        %% =============================================

        PD_roi = ...
            PDWF(y_start:y_end,...
               x_start:x_end,:);

        bfSIG_rf_roi = ...
            bfSIG_rf(y_start:y_end,...
                     x_start:x_end,:,:);

        fSIG_roi = ...
            fSIG(y_start:y_end,...
                 x_start:x_end,:,:);

        fSIG_rf_roi = ...
            fSIG_rf(y_start:y_end,...
                 x_start:x_end,:,:);

        %% =============================================
        % FEATURE EXTRACTION
        %% =============================================

        [nx_roi,ny_roi,npulse,nframe] = ...
            size(fSIG_roi);

        x1_map = zeros(nx_roi,ny_roi,nframe);

        x2_map = zeros(nx_roi,ny_roi,nframe);

        x1_frame = zeros(nframe,1);

        x2_frame = zeros(nframe,1);

        for n = 1:nframe

            fprintf('Feature frame %d / %d\n', ...
                n,nframe);

            for ix = 1:nx_roi

                for iy = 1:ny_roi

                    s = squeeze(fSIG_roi(ix,iy,:,n));

                    I = real(s);

                    Q = imag(s);

                    num_im = ...
                        sum( ...
                        I(2:end).*Q(1:end-1) - ...
                        I(1:end-1).*Q(2:end));

                    denom = ...
                        sum(I.^2 + Q.^2) + eps;

                    x1 = abs(num_im/denom);

                    num_re = ...
                        sum( ...
                        I(2:end).*I(1:end-1) + ...
                        Q(2:end).*Q(1:end-1));

                    R1_R0 = num_re/denom;

                    x2 = 1 - R1_R0;

                    x1_map(ix,iy,n) = x1;

                    x2_map(ix,iy,n) = x2;

                end
            end

            x1_frame(n) = ...
                mean(x1_map(:,:,n),'all');

            x2_frame(n) = ...
                mean(x2_map(:,:,n),'all');

        end

        %% =============================================
        % SAVE TWINKLING DATA
        %% =============================================

        save(fullfile(outputFolder,...
              [resultBase 'Metal_100Hz_roi.mat']),...
            'x1_map',...
            'x2_map',...
            'x1_frame',...
            'x2_frame',...
            '-v7.3');

        %% =============================================
        % SAVE ROI DATA
        %% =============================================

        ROI_Data = struct;

        ROI_Data.PD_roi = PD_roi;

        ROI_Data.bfSIG_rf_roi = bfSIG_rf_roi;

        ROI_Data.fSIG_roi = fSIG_roi;
        ROI_Data.fSIG_rf_roi=fSIG_rf_roi;

        ROI_Data.x_start = x_start;
        ROI_Data.x_end   = x_end;

        ROI_Data.y_start = y_start;
        ROI_Data.y_end   = y_end;

        save(fullfile(outputFolder,...
            [resultBase '_ROI_Data.mat']),...
            'ROI_Data','-v7.3');

        fprintf('\nSaved:\n%s\n',resultBase);

    end
end

fprintf('\n=====================================\n');
fprintf('ALL PROCESSING FINISHED\n');
fprintf('=====================================\n');

toc;
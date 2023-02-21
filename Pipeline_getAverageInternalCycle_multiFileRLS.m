%Pipeline_getAverageInternalCycle_multiFileRLS - examplary pipeline for processing
%rls file into contrast of a single period of activity of interest. E.g. it
%allows converting high-framerate data into a single cardiac cycle. Can also
%be used for vasomotion and other types of endogeneous perdiodic activity
%which is strong enough to be identified in the averaged signal. Can
%process multiple files in a sequence

% Authors: Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 17-Jan-2023

%------------- BEGIN CODE --------------
%% Pipeline configuration
%Clear the workspace to avoid potential conflicts.
clear
%Add all LSCI library folders to the matlab path
lsciLibraryFolder = 'C:\Users\AU707705\Dropbox\Work\GitHub\DDPLab-private\Light scattering';
addpath(genpath(lsciLibraryFolder));

settings.decimationSpace=4; %spatial decimation used to conserve memory in the pre-processing steps
settings.framesToAverage=1; %allows averaging multiple raw frames to artificially increase expsoure time
settings.contrastKernelS=5; %contrast kernel for spatial (sLSCI) processing method
settings.contrastKernelT=25; %contrast kernel for temporal (tLSCI) and lossless (ltLSCI) processing methods
settings.contrastKernelPreproc=5; %contrast kernel used in preprocessing (spatial)
settings.trustLimitsK=[0.0001,0.5];
settings.maxFrqIni=20; % initial max frequency of the activity of interest, Hz
settings.minFrqIni=3; % initial min frequency of the activity of interest, Hz
settings.rangeFrq=1/2; % relative frq range around the central frequency, Hz

settings.interpFactor=10; %Sets the number of points that will replace two consequitive values during the interpolation sequence.
settings.smoothCoef1=1/3; %in respect to minimum points per cycle value
settings.minPromCoef=1/2; % in respect to the std of the signal
settings.method={'sLSCIMM'};%,'tLSCIMM','ltLSCIMM'}; %Typically 'sLSCIMM' is recommended. For high quality data 'ltLSCIMM' will produce better results. Other options are 'tLSCIMM' and 'sLSCIMMM'.
% method refers to spatial, temporal or lossless contrast calculation,
% while the MM or MMM refers to minimum to minimum stretching or minimum to
% maximum + maximum to minimum stretching.

settings.enableRejectionModification=false; %allows manual correction of accepted pulses - disabled by default
settings.coeffsSTD=[2,2,2,2,2,2,2,2,2]; %pulses rejection coefficients relative to the feature standard deviation
settings.coeffsRel=0.25; %pulses rejection coefficients relative to the feature value
settings.coeffsAbs=5; %pulses rejection coefficients relative to the absolute feature value
settings.excludeFirstCycle=1; %always reject the first cycle


%Path to rls files (change indexes correspondingly to the number of files)
rawFileNames{1}='O:\HE_BFI-DATA\Mia\Pulsatility\2023111 - PSO01\2023111 - PSO01_awake.rls';
rawFileNames{2}='O:\HE_BFI-DATA\Mia\Pulsatility\2023111 - PSO01\2023111 - PSO01_iso.rls';
rawFileNames{3}='O:\HE_BFI-DATA\Mia\Pulsatility\2023111 - PSO01\2023111 - PSO01_k.rls';
rawFileNames{4}='O:\HE_BFI-DATA\Mia\Pulsatility\2023112 - PSM01\2023112 - PSM01_a.rls';
rawFileNames{5}='O:\HE_BFI-DATA\Mia\Pulsatility\2023112 - PSM01\2023112 - PSM01_i.rls';
rawFileNames{6}='O:\HE_BFI-DATA\Mia\Pulsatility\2023112 - PSM01\2023112 - PSM01_k.rls';
rawFileNames{7}='O:\HE_BFI-DATA\Mia\Pulsatility\2023112 - PSM02\2023112 - PSM02_a.rls';
rawFileNames{8}='O:\HE_BFI-DATA\Mia\Pulsatility\2023112 - PSM02\2023112 - PSM02_i.rls';
rawFileNames{9}='O:\HE_BFI-DATA\Mia\Pulsatility\2023112 - PSM02\2023112 - PSM02_k.rls';
rawFileNames{10}='O:\HE_BFI-DATA\Mia\Pulsatility\2023116 - PSY01\2023116 - PSY01_a_1.rls';
rawFileNames{11}='O:\HE_BFI-DATA\Mia\Pulsatility\2023116 - PSY01\2023116 - PSY01_i_1.2.rls';
rawFileNames{12}='O:\HE_BFI-DATA\Mia\Pulsatility\2023116 - PSY01\2023116 - PSY01_k_1.rls';
rawFileNames{13}='O:\HE_BFI-DATA\Mia\Pulsatility\2023117 - PSO04\2023117 - PSO04_a.rls';
rawFileNames{14}='O:\HE_BFI-DATA\Mia\Pulsatility\2023117 - PSO04\2023117 - PSO04_i.rls';
rawFileNames{15}='O:\HE_BFI-DATA\Mia\Pulsatility\2023117 - PSO04\2023117 - PSO04_k.rls';
rawFileNames{16}='O:\HE_BFI-DATA\Mia\Pulsatility\2023118 - PSY04\2023118 PSY04_a2.rls';
rawFileNames{17}='O:\HE_BFI-DATA\Mia\Pulsatility\2023118 - PSY04\2023118 PSY04_i.rls';
rawFileNames{18}='O:\HE_BFI-DATA\Mia\Pulsatility\2023118 - PSY04\2023118 PSY04_k.rls';
rawFileNames{19}='O:\HE_BFI-DATA\Mia\Pulsatility\2023118 - PSY05\2023118 - PSY05_a2.rls';
rawFileNames{20}='O:\HE_BFI-DATA\Mia\Pulsatility\2023118 - PSY05\2023118 - PSY05_i.rls';
rawFileNames{21}='O:\HE_BFI-DATA\Mia\Pulsatility\2023118 - PSY05\2023118 - PSY05_k.rls';
rawFileNames{22}='O:\HE_BFI-DATA\Mia\Pulsatility\2023119 - PSM05\2023119 - PSM05_a.rls';
rawFileNames{23}='O:\HE_BFI-DATA\Mia\Pulsatility\2023119 - PSM05\2023119 - PSM05_i.rls';
rawFileNames{24}='O:\HE_BFI-DATA\Mia\Pulsatility\2023119 - PSM05\2023119 - PSM05_k.rls';
rawFileNames{25}='O:\HE_BFI-DATA\Mia\Pulsatility\2023119 - PSM03\2023119 - PSM03_a.rls';
rawFileNames{26}='O:\HE_BFI-DATA\Mia\Pulsatility\2023119 - PSM03\2023119 - PSM03_a2.rls';
rawFileNames{27}='O:\HE_BFI-DATA\Mia\Pulsatility\2023119 - PSM03\2023119 - PSM03_i.rls';
rawFileNames{28}='O:\HE_BFI-DATA\Mia\Pulsatility\2023119 - PSM03\2023119 - PSM03_k.rls';
rawFileNames{29}='O:\HE_BFI-DATA\Mia\Pulsatility\2023118 - PSY05\2023118 - PSY05_a.rls';

for fidx=1:1:length(rawFileNames)
    tic
    clearvars results;
    close all
    settings.fName=rawFileNames{fidx};

    %read the raw file meta data
    fid = fopen(settings.fName, 'r');
    fseek(fid,0*1024,-1 );
    settings.sizeX=fread(fid,1,'*uint64');
    settings.sizeY=fread(fid,1,'*uint64');
    settings.sizeT=single(fread(fid,1,'*uint64'))-1;
    settings.exposureTime=single(fread(fid,1,'*uint64'));
    settings.fps=1000./settings.exposureTime; %converting time between frames to fps. DOES NOT PROVIDE AN ACCURATE FPS VALUES, RE-EVALUATED BELOW.
    settings.dataSize=1; %set to 1 by default for uint8 data type.
    settings.version=fread(fid,4,'*ubit8')';
    if strcmp(settings.version,'Ver.')
        settings.nVer = fread(fid,1,'*uint64');
        if nVer>1
            settings.dataSize=fread(fid,1,'*uint64');
        end
    end
    switch settings.dataSize
        case 1
            settings.dataType='*uint8';
        case 2
            settings.dataType='*uint16';
        otherwise
            error('Unindentified data type')
    end

    settings.sizeT=settings.sizeT-settings.framesToAverage+1;
    data=zeros(floor(settings.sizeY/settings.decimationSpace),floor(settings.sizeX/settings.decimationSpace),settings.sizeT,'single');
    meanI=zeros(1,settings.sizeT,'single');
    timeStamps=zeros(settings.sizeT,1,'int64');

    %move to the first timeStamp/frame location
    firstByte=30*1024;
    fseek(fid,firstByte,-1 );

    %perform data pre-processing, store spacially decimated information
    kernel=gpuArray(ones(settings.contrastKernelPreproc,settings.contrastKernelPreproc,1,'single'));
    for i=1:1:size(data,3)
        timeStamps(i)=fread(fid,1,'*uint64');
        tmp=gpuArray(fread(fid,settings.sizeX*settings.sizeY*settings.framesToAverage,settings.dataType));
        meanI(i)=mean(tmp(:));
        tmp=reshape(tmp,settings.sizeY,settings.sizeX,settings.framesToAverage);
        tmp=mean(single(tmp),3);
        tmp=stdfilt(tmp,kernel)./ convn(tmp,kernel,'same') * length(kernel(:));
        tmp=convn(tmp,ones(settings.decimationSpace),'same')./(settings.decimationSpace.*settings.decimationSpace);
        data(:,:,i)=gather(tmp(1:settings.decimationSpace:end,1:settings.decimationSpace:end));
    end




    imgK=mean(data,3);
    results.mask=imgK<settings.trustLimitsK(2) & imgK>settings.trustLimitsK(1) & ~isnan(imgK) & isfinite(imgK);
    results.mask=imerode(results.mask,[0,1,0;1,1,1;0,1,0]);
    data(isnan(data))=0;
    data(~isfinite(data))=0;
    data=reshape(data,size(data,1)*size(data,2),settings.sizeT);
    tsK=sum(data.*results.mask(:),1)./sum(results.mask(:));
    tsBFI=1./(tsK.*tsK);


    %correct the framerate
    settings.fps=1000./(mean(double(timeStamps(2:end)-timeStamps(1:end-1))));
    time=((1:1:length(tsK))-1)./settings.fps;

    %identify central frequency
    [fftPow,~,f]=getFFT(tsBFI,settings.fps,2.^(nextpow2(settings.fps/settings.maxFrqIni*50)),'cpu');
    [~,idx]=max(fftPow(:).*(f(:)>settings.minFrqIni).*(f(:)<settings.maxFrqIni));
    settings.centralFrq=f(idx);
    settings.minFrq=max(settings.centralFrq*(1-settings.rangeFrq),settings.minFrqIni);
    settings.maxFrq=min(settings.centralFrq*(1+settings.rangeFrq),settings.maxFrqIni);
    settings.meanPointsPerCycle=floor(settings.fps/settings.centralFrq);


    tsBFIF=smooth(tsBFI,max(floor(settings.meanPointsPerCycle.*settings.smoothCoef1),1),'loess');
    figure
    subplot(1,2,1)
    plot(time,tsBFIF)
    hold on
    plot(time,1./(tsK.*tsK))
    hold off
    settings.minProm=settings.minPromCoef.*std(tsBFIF);

    [~,locsMin]=findpeaks(-tsBFIF,'MinPeakDistance',floor(settings.fps/settings.maxFrq),'MinPeakProminence',settings.minProm);

    locsMax=zeros(length(locsMin)-1,1);
    for i=1:1:length(locsMin)-1
        [~,idx]=max(tsBFIF(locsMin(i):locsMin(i+1)));
        locsMax(i)=locsMin(i)+idx-1;

    end
    results.pulsesList=zeros(length(locsMax),3);
    for i=1:1:length(locsMin)-1
        results.pulsesList(i,:)=[locsMin(i),locsMax(i),locsMin(i+1)];
    end





    results.pulsesFeatures=zeros(size(results.pulsesList,1),9);
    for i=1:1:size(results.pulsesFeatures,1)
        results.pulsesFeatures(i,1)=mean(tsBFIF(squeeze(results.pulsesList(i,:))));
        results.pulsesFeatures(i,2)=std(tsBFIF(squeeze(results.pulsesList(i,:))));
        results.pulsesFeatures(i,3)=max(tsBFIF(squeeze(results.pulsesList(i,:))))-min(tsBFIF(squeeze(results.pulsesList(i,:))));
        results.pulsesFeatures(i,4)=max(squeeze(results.pulsesList(i,:)))-min(squeeze(results.pulsesList(i,:)));
        results.pulsesFeatures(i,5)=abs(tsBFIF(squeeze(results.pulsesList(i,1)))-tsBFIF(squeeze(results.pulsesList(i,3))));
        results.pulsesFeatures(i,6)=sum((tsBFIF(results.pulsesList(i,1)+1:results.pulsesList(i,2))-tsBFIF(results.pulsesList(i,1):results.pulsesList(i,2)-1))<0);
        results.pulsesFeatures(i,7)=sum((tsBFIF(results.pulsesList(i,2)+1:results.pulsesList(i,3))-tsBFIF(results.pulsesList(i,2):results.pulsesList(i,3)-1))>0);
        results.pulsesFeatures(i,8)=mean((tsBFIF(results.pulsesList(i,2)+1:results.pulsesList(i,3))-tsBFIF(results.pulsesList(i,2):results.pulsesList(i,3)-1)));
        results.pulsesFeatures(i,9)=mean(meanI(results.pulsesList(i,1)+1:results.pulsesList(i,3)));
    end


    results.pulsesToReject=zeros(11,size(results.pulsesList,1));
    for i=1:1:size(results.pulsesList,1)

        for ii=1:1:size(results.pulsesFeatures,2)
            if abs(results.pulsesFeatures(i,ii)-median(squeeze(results.pulsesFeatures([1:i-1,i+1:end],ii))))>settings.coeffsSTD(ii)*std(squeeze(results.pulsesFeatures([1:i-1,i+1:end],ii)))
                results.pulsesToReject(ii,i)=1;
            end
        end

        if results.pulsesFeatures(i,5)/results.pulsesFeatures(i,3)>settings.coeffsRel(1)
            results.pulsesToReject(10,i)=1;
        end
        if (settings.fps/results.pulsesFeatures(i,4))<settings.minFrq || (settings.fps/results.pulsesFeatures(i,4))>settings.maxFrq
            results.pulsesToReject(11,i)=1;
        end

        if any(strcmp(settings.method,'tLSCIMM'))
            if results.pulsesList(i,3)>(settings.sizeT-settings.contrastKernelT)
                results.pulsesToReject(12,i)=1;
            end
        end

    end

    if settings.excludeFirstCycle==1
        results.pulsesToReject(1:end,1)=1;
    end


    %reject/accept based on number of consequtive accepts/rejects
    a=max(results.pulsesToReject,[],1);
    a=movmean(a,[settings.coeffsAbs(1),settings.coeffsAbs(1)]);
    results.pulsesToReject(13,:)=round(a);
    %b=diff([0 find(diff(a)) numel(a)]);
    %c = arrayfun(@(x) ones(1,x).*x, b , 'un',0);
    %c=cat(2,c{:});
    %results.pulsesToReject(13,:)=c<settings.coeffsAbs(1);


    subplot(1,2,2)
    hold on
    for i=1:1:size(results.pulsesList,1)
        plot([time(results.pulsesList(i,1)),time(results.pulsesList(i,3))],[tsBFIF(results.pulsesList(i,1)),tsBFIF(results.pulsesList(i,3))],'ok')
        if sum(results.pulsesToReject(:,i))==0
            plot(time(results.pulsesList(i,1):results.pulsesList(i,3)),tsBFIF(results.pulsesList(i,1):results.pulsesList(i,3)),'g')
        else
            plot(time(results.pulsesList(i,1):results.pulsesList(i,3)),tsBFIF(results.pulsesList(i,1):results.pulsesList(i,3)),'r')
        end
        for ii=1:1:size(results.pulsesToReject,1)
            if results.pulsesToReject(ii,i)==0
                plot([time(results.pulsesList(i,1)),time(results.pulsesList(i,3))],[min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(results.pulsesToReject,1)/2,min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(results.pulsesToReject,1)/2],'g');
            else
                plot([time(results.pulsesList(i,1)),time(results.pulsesList(i,3))],[min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(results.pulsesToReject,1)/2,min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(results.pulsesToReject,1)/2],'r');
            end

        end
    end
    hold off
    if settings.enableRejectionModification
        [x,~]=ginput();
        for i=1:1:length(x)
            results.pulsesToReject(:,x(i)>=time(results.pulsesList(:,1)) & x(i)<time(results.pulsesList(:,3)))=1-max(results.pulsesToReject(:,x(i)>=time(results.pulsesList(:,1)) & x(i)<time(results.pulsesList(:,3))));
        end

        subplot(1,2,2)
        hold on
        for i=1:1:size(results.pulsesList,1)
            plot([time(results.pulsesList(i,1)),time(results.pulsesList(i,3))],[tsBFIF(results.pulsesList(i,1)),tsBFIF(results.pulsesList(i,3))],'ok')
            if sum(results.pulsesToReject(:,i))==0
                plot(time(results.pulsesList(i,1):results.pulsesList(i,3)),tsBFIF(results.pulsesList(i,1):results.pulsesList(i,3)),'g')
            else
                plot(time(results.pulsesList(i,1):results.pulsesList(i,3)),tsBFIF(results.pulsesList(i,1):results.pulsesList(i,3)),'r')
            end
            for ii=1:1:size(results.pulsesToReject,1)
                if results.pulsesToReject(ii,i)==0
                    plot([time(results.pulsesList(i,1)),time(results.pulsesList(i,3))],[min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(results.pulsesToReject,1)/2,min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(results.pulsesToReject,1)/2],'g');
                else
                    plot([time(results.pulsesList(i,1)),time(results.pulsesList(i,3))],[min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(results.pulsesToReject,1)/2,min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(results.pulsesToReject,1)/2],'r');
                end

            end
        end
        hold off
    end
    xlabel('Time,s')
    ylabel('BFI')
    drawnow
    toc
    ax = gca;
    exportgraphics(ax,[settings.fName,'1.jpg'],'Resolution',600);

    results.pulsesListFinal=results.pulsesList;
    for i=size(results.pulsesList,1):-1:1
        if sum(results.pulsesToReject(:,i))>0
            results.pulsesListFinal(i,:)=[];
        end
    end
    clearvars data;

    results.descendTimePts=round(median(results.pulsesListFinal(:,3)-results.pulsesListFinal(:,2)));
    results.ascendTimePts=round(median(results.pulsesListFinal(:,2)-results.pulsesListFinal(:,1)));
    results.cycleTimePts=results.ascendTimePts+results.descendTimePts-1;


    %SLSCI method Min Max Min
    if any(strcmp(settings.method,'sLSCIMMM'))
        dataASLSCI=zeros(results.ascendTimePts*settings.interpFactor,settings.sizeY,settings.sizeX,'single');
        dataDSLSCI=zeros(results.descendTimePts*settings.interpFactor,settings.sizeY,settings.sizeX,'single');

        kernel=ones(settings.contrastKernelS,settings.contrastKernelS,1);
        for i=1:1:size(results.pulsesListFinal,1)
            data=readRLS(settings.fName,results.pulsesListFinal(i,1)-1,results.pulsesListFinal(i,3)-results.pulsesListFinal(i,1)+settings.framesToAverage);
            data=movmean(data,[0,settings.framesToAverage-1],3,'Endpoints','discard');
            data= stdfilt (data,kernel)./ convn(data,kernel,'same') * length(kernel(:));

            ascendIdxs=1:(results.pulsesListFinal(i,2)-results.pulsesListFinal(i,1)+1);
            descendIdxs=(results.pulsesListFinal(i,2)-results.pulsesListFinal(i,1)+1):(results.pulsesListFinal(i,3)-results.pulsesListFinal(i,1)+1);

            data=permute(data,[3,1,2]);
            dataASLSCI=dataASLSCI+interp1(data(ascendIdxs,:,:),linspace(1,length(ascendIdxs),results.ascendTimePts*settings.interpFactor));
            dataDSLSCI=dataDSLSCI+interp1(data(descendIdxs,:,:),linspace(1,length(descendIdxs),results.descendTimePts*settings.interpFactor));



        end
        dataASLSCI=permute(dataASLSCI,[2,3,1]);
        dataDSLSCI=permute(dataDSLSCI,[2,3,1]);
        results.data=cat(3,dataASLSCI(:,:,1:end-1),dataDSLSCI)./size(results.pulsesListFinal,1);
        results.time=(0:1:size(dataCycle,3)-1)./settings.fps./settings.interpFactor;
        dataALSCI=[];
        dataDLSCI=[];
        tmp=strsplit(settings.fName,'\');
        h=figure;
        subplot(2,1,1)
        img=squeeze(mean(results.data,3));
        imagesc(img);
        caxis([prctile(img(:),5),prctile(img(:),99)]);
        axis image
        subplot(2,1,2)
        plot(results.time,1./(squeeze(mean(results.data,[1,2]))).^2);
        xlabel('Time,s')
        ylabel('BFI')
        sgtitle(['sLSCIMMM ',strrep(tmp{end},'_',' ')]);
        drawnow
        disp('Saving the results sLSCIMMM');
        tmp=strrep(settings.fName,'.rls','_cycle_sLSCIMMM.mat');
        if contains(tmp,'__cycle_sLSCIMMM.mat')
            save(tmp,'results','settings','-v7.3');
            disp('Saving complete');
        else
            disp('Saving aborted, wrong file name')
        end
        toc
        clearvars data dataASLSCI dataDSLSCI;

    end
    if any(strcmp(settings.method,'sLSCIMM'))
        %SLSCI method Min Min
        dataCycle=zeros(results.cycleTimePts*settings.interpFactor,settings.sizeY,settings.sizeX,'single');
        kernel=ones(settings.contrastKernelS,settings.contrastKernelS,1);
        for i=1:1:size(results.pulsesListFinal,1)
            data=readRLS(settings.fName,results.pulsesListFinal(i,1)-1,results.pulsesListFinal(i,3)-results.pulsesListFinal(i,1)+settings.framesToAverage);
            data=movmean(data,[0,settings.framesToAverage-1],3,'Endpoints','discard');
            data= stdfilt (data,kernel)./ convn(data,kernel,'same') * length(kernel(:));
            data=permute(data,[3,1,2]);
            idxs=1:1:(results.pulsesListFinal(i,3)-results.pulsesListFinal(i,1)+1);
            dataCycle=dataCycle+interp1(data(idxs,:,:),linspace(1,length(idxs),results.cycleTimePts*settings.interpFactor));

        end
        dataCycle=permute(dataCycle,[2,3,1]);
        results.data=dataCycle./size(results.pulsesListFinal,1);
        results.time=(0:1:size(dataCycle,3)-1)./settings.fps./settings.interpFactor;
        tmp=strsplit(settings.fName,'\');
        h=figure;
        subplot(2,1,1)
        img=squeeze(mean(results.data,3));
        imagesc(img);
        caxis([prctile(img(:),5),prctile(img(:),99)]);
        axis image
        subplot(2,1,2)
        plot(results.time,1./(squeeze(mean(results.data,[1,2]))).^2);
        xlabel('Time,s')
        ylabel('BFI')
        sgtitle(['sLSCIMM ',strrep(tmp{end},'_',' ')]);
        drawnow
        disp('Saving the results sLSCIMM');
        tmp=strrep(settings.fName,'.rls','_cycle_sLSCIMM.mat');
        if contains(tmp,'_cycle_sLSCIMM.mat')
            save(tmp,'results','settings','-v7.3');
            disp('Saving complete');
        else
            disp('Saving aborted, wrong file name')
        end
        toc
        clearvars data dataCycle;
    end

    if any(strcmp(settings.method,'tLSCIMM'))
        %TLSCI method Min Min
        %use with care TLSCI might be contaminated by the rejected pulses
        dataCycle=zeros(results.cycleTimePts*settings.interpFactor,settings.sizeY,settings.sizeX,'single');
        for i=1:1:size(results.pulsesListFinal,1)
            data=readRLS(settings.fName,results.pulsesListFinal(i,1)-1-(floor(settings.contrastKernelT/2*settings.framesToAverage)),results.pulsesListFinal(i,3)-results.pulsesListFinal(i,1)+settings.framesToAverage+(settings.contrastKernelT*settings.framesToAverage-1));
            data=movmean(data,[0,settings.framesToAverage-1],3,'Endpoints','discard');
            data= movstd (data,[0,settings.contrastKernelT*settings.framesToAverage-1],0,3,'Endpoints','discard')...
                ./ movmean(data,[0,settings.contrastKernelT*settings.framesToAverage-1],3,'Endpoints','discard');
            data=permute(data,[3,1,2]);
            idxs=1:1:(results.pulsesListFinal(i,3)-results.pulsesListFinal(i,1)+1);
            dataCycle=dataCycle+interp1(data(idxs,:,:),linspace(1,length(idxs),results.cycleTimePts*settings.interpFactor));

        end
        dataCycle=permute(dataCycle,[2,3,1]);
        results.data=dataCycle./size(results.pulsesListFinal,1);
        results.time=(0:1:size(dataCycle,3)-1)./settings.fps./settings.interpFactor;
        tmp=strsplit(settings.fName,'\');
        h=figure;
        subplot(2,1,1)
        img=squeeze(mean(results.data,3));
        imagesc(img);
        caxis([prctile(img(:),5),prctile(img(:),99)]);
        axis image
        subplot(2,1,2)
        plot(results.time,1./(squeeze(mean(results.data,[1,2]))).^2);
        xlabel('Time,s')
        ylabel('BFI')
        sgtitle(['tLSCIMM ',strrep(tmp{end},'_',' ')]);
        drawnow
        disp('Saving the results tLSCIMM');
        tmp=strrep(settings.fName,'.rls','_cycle_tLSCIMM.mat');
        if contains(tmp,'_cycle_tLSCIMM.mat')
            save(tmp,'results','settings','-v7.3');
            disp('Saving complete');
        else
            disp('Saving aborted, wrong file name')
        end
        toc
        clearvars data dataCycle;
    end

    if any(strcmp(settings.method,'ltLSCIMM'))
        %LTLSCI method Min Min
        dataCycle=zeros(results.cycleTimePts*settings.interpFactor,settings.sizeY,settings.sizeX,'single');
        cyclesDurList=(results.pulsesListFinal(:,3)-results.pulsesListFinal(:,1)+1);
        cyclesDur=unique(cyclesDurList);
        cyclesN=0;

        for i=1:1:length(cyclesDur)
            pulsesIdxs=find(cyclesDurList==cyclesDur(i));
            if length(pulsesIdxs)>=settings.contrastKernelT
                data=zeros(settings.sizeY,settings.sizeX,cyclesDur(i),length(pulsesIdxs),'single');
                for ii=1:1:length(pulsesIdxs)
                    data(:,:,:,ii)=movmean(readRLS(settings.fName,results.pulsesListFinal(pulsesIdxs(ii),1)-1,cyclesDur(i)+settings.framesToAverage-1),[0,settings.framesToAverage-1],3,'Endpoints','discard');
                end
                data=squeeze(std(data,0,4)./mean(data,4));
                data=permute(data,[3,1,2]);
                dataCycle=dataCycle+length(pulsesIdxs).*interp1(data,linspace(1,cyclesDur(i),results.cycleTimePts*settings.interpFactor));
                cyclesN=cyclesN+length(pulsesIdxs);

            end
        end
        dataCycle=permute(dataCycle,[2,3,1]);
        results.data=dataCycle./cyclesN;
        results.time=(0:1:size(dataCycle,3)-1)./settings.fps./settings.interpFactor;
        tmp=strsplit(settings.fName,'\');
        h=figure;
        subplot(2,1,1)
        img=squeeze(mean(results.data,3));
        imagesc(img);
        caxis([prctile(img(:),5),prctile(img(:),99)]);
        axis image
        subplot(2,1,2)
        plot(results.time,1./(squeeze(mean(results.data,[1,2]))).^2);
        xlabel('Time,s')
        ylabel('BFI')
        sgtitle(['ltLSCIMM ',strrep(tmp{end},'_',' ')]);
        drawnow
        disp('Saving the results ltLSCIMM');
        tmp=strrep(settings.fName,'.rls','_cycle_ltLSCIMM.mat');
        if contains(tmp,'_cycle_ltLSCIMM.mat')
            save(tmp,'results','settings','-v7.3');
            disp('Saving complete');
        else
            disp('Saving aborted, wrong file name')
        end
        toc
        clearvars data dataCycle;

    end
    exportgraphics(h,[settings.fName,'2.jpg'],'Resolution',600);

end
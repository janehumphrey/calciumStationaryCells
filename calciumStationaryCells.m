%%
% NAME: CALCIUM (STATIONARY CELLS)
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

tic;
close all;
clear variables;
clc;

% USER INPUT.
inputDir = '';  % Folder containing a single .tif file.
bgFile = '';  % Background file (optional).
resultsDir = '';  % Results folder.
Info.cellDiam = 20;  % Cell diameter, in um.
Info.pixelSize = 16/20/1.5;  % Pixel size, in um.
Info.interval = 1;  % Time between frames, in s.
Info.minInt = 0.05;  % Intensity threshold for cell selection (SDs above mean). After running check Peaks.tif to adjust.
Info.smoothing = 5;  % Standard deviation for smoothing of traces, in s.
Info.minGrad = 0.01;  % Threshold (intensity gradient) defining a calcium spike, in /s. After running check intensity
% traces to adjust.
Info.minHeight = 1.25;  % Minimum intensity of a calcium spike, relative to baseline (optional).
Info.minWidth = 0;  % Minimum duration of a calcium spike, in s (optional).
%

functionDir = [pwd,filesep,'Functions'];
if isfolder(functionDir)
    addpath(genpath(functionDir));
    addpath(genpath([pwd,filesep,'Extra_functions']));
else
    addpath(genpath([fileparts(pwd),filesep,'Functions']));
    addpath(genpath([fileparts(pwd),filesep,'Extra_functions']));
end
if ~isfolder(inputDir)
    error('Directory doesn''t exist.');
end

[fileList,fileExt,nFiles] = selectFiles(inputDir,'tif');
if nFiles==1
    file = [inputDir,filesep,fileList{1},fileExt{1}];
else
    error('Too many files in directory.');
end

resultsDir = createDir(resultsDir,inputDir,'Calcium');
trueDir = [resultsDir,filesep,'Calcium_release'];
mkdir(trueDir);
falseDir = [resultsDir,filesep,'No_calcium_release'];
mkdir(falseDir);

radiusPix = round(Info.cellDiam/Info.pixelSize/2);
[fovRaw,width,height,nFrames] = readStack(file);
fovDp = double(fovRaw);
if isempty(bgFile)
    bgMean = min(fovDp(:));
else
    bgRaw = readStack(bgFile);
    bgMean = mean(bgRaw,3);
end
fovNoBg = fovDp-bgMean;

fovAvg = mean(fovNoBg,3);
fovFilt = imgaussfilt(fovAvg,radiusPix*10);
missing2d = fovFilt<=1;
missing3d = repmat(missing2d,[1,1,nFrames]);
fovFlat = fovNoBg./fovFilt;
fovFlat(missing3d) = 0;
gaussSd = round(radiusPix/2);
fovSmooth = imgaussfilt(fovFlat,gaussSd);

fovMean = mean(fovSmooth(:));
fovSd = mean(fovSmooth(:));
peakThresh = fovMean+Info.minInt*fovSd;
frameNos = (1:nFrames)';
time = frameNos*Info.interval;

coord = cell(nFrames,1);
peaks = cell(nFrames,1);
for iFrame = 1:nFrames
    coord{iFrame} = pkfnd(fovSmooth(:,:,iFrame),peakThresh,radiusPix*1.5);
    nPeaks = size(coord{iFrame},1);
    peaks{iFrame} = NaN(nPeaks,3);
    peaks{iFrame}(:,1:2) = coord{iFrame};
    peaks{iFrame}(:,3) = iFrame;
end
peaks = vertcat(peaks{:});

trackingParam.mem = 5;
trackingParam.dim = 2;
trackingParam.good = nFrames*0.51;
trackingParam.quiet = 1;
trackMat = track(peaks,radiusPix,trackingParam);

nTracks = trackMat(end,4);
goodTracks = true(nTracks,1);
Results.xPos = NaN(nTracks,1);
Results.yPos = NaN(nTracks,1);
for iTrack = 1:nTracks
    trackValues = trackMat(trackMat(:,4)==iTrack,:);
    xValues = trackValues(:,1);
    yValues = trackValues(:,2);
    Results.xPos(iTrack) = round(nanmean(xValues));     
    Results.yPos(iTrack) = round(nanmean(yValues));
    xSd = nanstd(xValues);
    ySd = nanstd(yValues);
    if xSd>radiusPix/2||ySd>radiusPix/2
        goodTracks(iTrack) = false;
    end
end
Results.xPos = Results.xPos(goodTracks);
Results.yPos = Results.yPos(goodTracks);
nCells = sum(goodTracks);
cellNos = (1:nCells)';

if nCells==0
    error('No cells selected.');
end

normTraces = NaN(nFrames,nCells);
gradTraces = NaN(nFrames,nCells);
spikeLoc = cell(nCells,1);
areas = cell(nCells,1);
cellId = cell(nCells,1);
Spikes.no = cell(nCells,1);
Spikes.time = cell(nCells,1);
Spikes.height = cell(nCells,1);
Spikes.width = cell(nCells,1);
Spikes.area = cell(nCells,1);
Results.caRelease = NaN(nCells,1);
Results.nSpikes = NaN(nCells,1);
Results.initialTime = NaN(nCells,1);
Results.maxHeight = NaN(nCells,1);
Results.maxWidth = NaN(nCells,1);
Results.combArea = NaN(nCells,1);

se = strel('disk',radiusPix,0);
warning('off','signal:findpeaks:largeMinPeakHeight');

for iCell = 1:nCells
    goodSpikes = [];
    roiX = Results.xPos(iCell)-radiusPix:Results.xPos(iCell)+radiusPix;
    roiY = Results.yPos(iCell)-radiusPix:Results.yPos(iCell)+radiusPix;
    goodWidth = roiX>0&roiX<=width;
    goodHeight = roiY>0&roiY<=height;
    roiBox = fovSmooth(roiY(goodHeight),roiX(goodWidth),:);
    roiDisk = roiBox.*se.Neighborhood(goodHeight,goodWidth);
    intTrace = squeeze(nanmean(nanmean(roiDisk,1),2));
    
    smoothTrace = smoothdata(intTrace,'gaussian',Info.smoothing*5/Info.interval);
    tempSmooth = smoothTrace;
    sortedTrace = sort(tempSmooth);
    cutoff = round(nFrames/4);
    refSmooth = sortedTrace(1:cutoff);
    baselineEst = mean(refSmooth);
    normTraces(:,iCell) = smoothTrace/baselineEst;
    gradTraces(:,iCell) = gradient(normTraces(:,iCell));
    tempGrad = gradTraces(:,iCell);
    
    [peakHeight,peakLoc,~,peakProm] = findpeaks(tempGrad,'minPeakHeight',Info.minGrad);
    spikeStart = peakLoc(peakProm>0.75*peakHeight);
    if ~isempty(spikeStart)
        heightAtStart = normTraces(spikeStart,iCell);
        nSpikes = size(spikeStart,1);
        overlapping = false(nSpikes,1);
        spikeEnd = NaN(nSpikes,1);
        spikeInt = NaN(nFrames,nSpikes);
        for iSpike = 1:nSpikes
            if max(spikeStart(iSpike)<spikeEnd(:))
                overlapping(iSpike) = true;
            end
            afterSpikeStart = frameNos>spikeStart(iSpike);
            belowStartHeight = normTraces(:,iCell)<heightAtStart(iSpike);
            if ~isempty(find(afterSpikeStart&belowStartHeight,1))
                spikeEnd(iSpike) = find(afterSpikeStart&belowStartHeight,1);
            else
                spikeEnd(iSpike) = nFrames;
            end
            spikeInt(spikeStart(iSpike):spikeEnd(iSpike),iSpike) = normTraces(spikeStart(iSpike):spikeEnd(iSpike),...
                iCell);
        end
        spikeHeight = max(spikeInt,[],1)';
        spikeWidth = (spikeEnd-spikeStart)*(Info.interval);
        goodHeight = spikeHeight>=Info.minHeight;
        goodWidth = spikeWidth>=Info.minWidth;
        goodSpikes = goodHeight&goodWidth&~overlapping;
    end
    if any(goodSpikes)==1
        spikeLoc{iCell} = spikeStart(goodSpikes);
        Spikes.time{iCell} = spikeStart(goodSpikes)*Info.interval;
        Spikes.height{iCell} = spikeHeight(goodSpikes);
        Spikes.width{iCell} = spikeWidth(goodSpikes);
        areas{iCell} = spikeInt(:,goodSpikes);
        Spikes.area{iCell} = nansum(areas{iCell}-1,1)';
        
        Results.caRelease(iCell) = 1;
        Results.nSpikes(iCell) = size(spikeLoc{iCell},1);
        Results.initialTime(iCell) = Spikes.time{iCell}(1);
        Results.maxHeight(iCell) = max(Spikes.height{iCell});
        Results.maxWidth(iCell) = max(Spikes.width{iCell});
        Results.combArea(iCell) = sum(Spikes.area{iCell});
        
        cellId{iCell} = ones(Results.nSpikes(iCell),1)*iCell;
        Spikes.no{iCell} = (1:Results.nSpikes(iCell))';
    else
        Results.caRelease(iCell) = 0;
    end
end

cellsTrue = cellNos(Results.caRelease==1);
cellsFalse = cellNos(Results.caRelease==0);
maxDigits = size(num2str(cellNos(end)),2);
format = ['%0.',num2str(maxDigits),'u'];
cellList = num2str(cellNos,format);

prepareFigure(width/height,width/height);
ax = prepareAxes('image',[],[],[1-10,width+10],[1-10,height+10]);
xCorners = [1-10,width+10,width+10,1-10]';
yCorners = [1-10,1-10,height+10,height+10]';
lastFrame = fovSmooth(:,:,end);
minInt = quantile(lastFrame(:),0.2);    
maxInt = quantile(lastFrame(:),0.99);    
lutLimits = [minInt,maxInt];
if radiusPix/height>0.01
    markerSize = 100;
    fontSize = 11;
else
    markerSize = 60;
    fontSize = 8;
end
peaksFile = [resultsDir,filesep,'Peaks.tif'];
pointsFile = [resultsDir,filesep,'Cells.tif'];
labelsFile = [resultsDir,filesep,'Labelled_cells.tif'];

for iFrame = 100:100:nFrames
    patch(xCorners,yCorners,[0,0.5,1],'EdgeColor','none');
    imagesc(fovSmooth(:,:,iFrame),lutLimits);
    plotData('points',coord{iFrame}(:,1),coord{iFrame}(:,2),[0,0.5,1],[],[],[],'none',markerSize);
    export_fig(peaksFile,'-tif','-append');
    ax.Children.delete;
end

for iFrame = 10:10:nFrames
    patch(xCorners,yCorners,[1,0.5,0],'EdgeColor','none');
    imagesc(fovSmooth(:,:,iFrame),lutLimits);

    xTrue = Results.xPos(Results.caRelease==1);
    yTrue = Results.yPos(Results.caRelease==1);
    xFalse = Results.xPos(Results.caRelease==0);
    yFalse = Results.yPos(Results.caRelease==0);
    trueStr = num2str(cellsTrue);
    falseStr = num2str(cellsFalse);
    trueArray = cellstr(trueStr);
    falseArray = cellstr(falseStr);
    
    plotData('points',xTrue,yTrue,[],[],[],[],'none',markerSize);
    plotData('points',xFalse,yFalse,[1,0,0],[],[],[],'none',markerSize);
    export_fig(pointsFile,'-tif','-append');
    ax.Children(1:end-2).delete;
    
    plotData('text',xTrue,yTrue,[],[],[],[],[],[],trueArray,fontSize);
    plotData('text',xFalse,yFalse,[1,0,0],[],[],[],[],[],falseArray,fontSize);
    export_fig(labelsFile,'-tif','-append');
    ax.Children.delete;
end
close;

binEdges = [-Inf;time];
ordered = sort(Results.initialTime,1);
binned = cumsum(histcounts(ordered,binEdges))';
triggerTime = binned/nCells*100;

prepareFigure;
prepareAxes('plot','time (s)','triggered cells (%)',[0,nFrames],[0,100]);
[~,index] = unique(triggerTime,'stable');
timeX = time([index;nFrames]);
timeY = triggerTime([index;nFrames]);
plotData('line',timeX,timeY);
timeBase = [resultsDir,filesep,'Triggering'];
saveImage(timeBase);
close;

prepareFigure(0.75);
subplot(2,1,1);
normMax = max(normTraces(:));
if normMax>nanmean(normTraces(:))*10
    normMax = nanmean(normTraces(:))*10;
end
topAx = prepareAxes('plot','time (s)','normalised intensity',[0,time(end)],[0,normMax]);
subplot(2,1,2);
gradMax = max(gradTraces(:));
if gradMax>Info.minGrad*10
    gradMax = Info.minGrad*10;
elseif gradMax<Info.minGrad*2
    gradMax = Info.minGrad*2;
end
gradLabel = ['intensity gradient (',getUnit('perS'),')'];
bottomAx = prepareAxes('plot','time (s)',gradLabel,[0,time(end)],[-gradMax,gradMax]);
gradThresh = [Info.minGrad,Info.minGrad];
plotData('line',[0,time(end)],gradThresh,[0,0,0],':',2);

for iCell = 1:nCells
    subplot(2,1,1);
    if Results.caRelease(iCell)==1
        for iSpike = 1:Results.nSpikes(iCell)
            patches = area(time,areas{iCell}(:,iSpike),1,'LineStyle',':','LineWidth',2,'FaceColor',[0.85,0.85,0.85]);
            patches.BaseLine.Visible = 'off';
        end
    end
    if Results.caRelease(iCell)==1
        plotData('line',time,normTraces(:,iCell));
    else
        plotData('line',time,normTraces(:,iCell),[1,0,0]);
    end
    subplot(2,1,2);
    if Results.caRelease(iCell)==1
        plotData('line',time,gradTraces(:,iCell));
    else
        plotData('line',time,gradTraces(:,iCell),[1,0,0]);
    end
    if Results.caRelease(iCell)==1
        spikeTime = time(spikeLoc{iCell});
        spikeGrad = gradTraces(spikeLoc{iCell},iCell);
        plotData('points',spikeTime,spikeGrad,[0,0,0],[],[],[],[],50);
        imageBase = [trueDir,filesep,'Cell_',cellList(iCell,:)];
    else
        imageBase = [falseDir,filesep,'Cell_',cellList(iCell,:)];
    end
    saveImage(imageBase);
    topAx.Children.delete;
    bottomAx.Children(1:end-1).delete;
end
close;

excel = actxserver('Excel.Application');

intBase = repmat('Cell ',nCells,1);
intMat = [intBase,num2str(cellNos)];
intTitles = cellstr(intMat)';
intValues = num2cell(normTraces);
intText = [{'Time (s)'},intTitles;num2cell(time),intValues];
intFile = [resultsDir,filesep,'Traces.xlsx'];
createSpreadsheet(excel,intFile,intText);

timeText = [{'Time','Triggered cells (%)'};num2cell(timeX),num2cell(timeY)];
timeFile = [resultsDir,filesep,'Triggering.xlsx'];
createSpreadsheet(excel,timeFile,timeText,0);

cellDiamStr = ['Cell diameter (',getUnit('um'),')'];
pixelSizeStr = ['Pixel size (',getUnit('um'),')'];
minGradStr = ['Min gradient (',getUnit('perS'),')'];
infoValues = struct2cell(Info)';
infoText = [{cellDiamStr,pixelSizeStr,'Interval (s)','Intensity threshold (SDs)','Smoothing (SDs)',minGradStr,...
    'Min spike height','Min spike duration (s)'};infoValues];
infoFile = [resultsDir,filesep,'Info.xlsx'];
createSpreadsheet(excel,infoFile,infoText);

cellId = vertcat(cellId{:});
spikesFields = fieldnames(Spikes);
for iField = 1:size(spikesFields,1)
    fieldName = spikesFields{iField};
    Spikes.(fieldName) = vertcat(Spikes.(fieldName){:});
    Spikes.(fieldName) = round(Spikes.(fieldName),3,'significant');
end
spikesValues = cell2mat(struct2cell(Spikes)');
spikesMean = round(nanmean(spikesValues,1),3,'significant');
spikesSd = round(nanstd(spikesValues,0,1),3,'significant');
if ~isempty(spikesValues)
    spikesHeadings = {'Cell','Spike','Time','Spike height','Spike duration (s)','Integrated intensity'};
    spikesText = [spikesHeadings;num2cell(cellId),num2cell(spikesValues);{'Mean'},num2cell(spikesMean);{'SD'},...
        num2cell(spikesSd)];
    spikesFile = [resultsDir,filesep,'Spikes.xlsx'];
    createSpreadsheet(excel,spikesFile,spikesText,2);
end

resultsFields = fieldnames(Results);
for iField = 1:size(resultsFields,1)
    fieldName = resultsFields{iField};
    Results.(fieldName) = round(Results.(fieldName),3,'significant');
end
resultsValues = cell2mat(struct2cell(Results)');
resultsMean = round(nanmean(resultsValues,1),3,'significant');
resultsSd = round(nanstd(resultsValues,0,1),3,'significant');
resultsHeadings = {'Cell','X position','Y position','Calcium release','No of spikes','Time of first spike',...
    'Max spike height','Max spike duration (s)','Total integrated intensity (all spikes)'};
resultsText = [resultsHeadings;num2cell(cellNos),num2cell(resultsValues);{'Mean'},num2cell(resultsMean);{'SD'},...
    num2cell(resultsSd)];
resultsFile = [resultsDir,filesep,'Results.xlsx'];
createSpreadsheet(excel,resultsFile,resultsText,2);

excel.Quit;
excel.delete;

toc;
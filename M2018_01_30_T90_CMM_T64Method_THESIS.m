% Preamble.
clc
clear
tic;
datetime('now')

%%

% General test information.

ExperimentNumber = 'T90';
ExperimentDateCode = '2018-01-30';
ExperimentTemperature = '24';

% Inputs.
localStructure = '\';
RootFile = ['\\ubccmps-sfpp01.ead.ubc.ca\Comp\Users\Stewart, Andrew\Output\Experiments\' ExperimentNumber '\' ExperimentDateCode '_' ExperimentNumber '_CMM\'];

%%

% General information on file structures.
Addendum = [ExperimentDateCode '_' ExperimentNumber '_CMM_' 'Point_Cloud_IndividualWrinkle_Output' localStructure];
RootFileOutput = [RootFile Addendum];
DataCSVFolder = [ExperimentDateCode '_' ExperimentNumber '_CMM_Point_Cloud'];
FileNameTerminator = '.xyz';

%%

% Start and endpoints of the test point cloud files to import.
StartRunFiles = 1;
EndRunFiles = 37;

%%

% Importing and cleaning of data.

d = dir([RootFile DataCSVFolder localStructure '*' FileNameTerminator]);
names = {d.name};

for i = 1:length(names)
    timeSignature{i} = datetime([names{i}(1:10) names{i}(end-8:end-4)],'InputFormat','yyyy-MM-dd_HHmm');
end

% Which file to set as the start test file for calculating the entire test duration.
startFrameTime = 1;

for i = 1:EndRunFiles-StartRunFiles+1
    RawImport{i} = [RootFile DataCSVFolder localStructure char(names(i-1+StartRunFiles))];
end

% Imports the run file.  Header can be removed if required.
for i = 1:EndRunFiles-StartRunFiles+1
    
    fileID{i} = fopen(RawImport{i},'r');
    dataArray{i} = textscan(fileID{i},'%f %f %f','Delimiter',' ','MultipleDelimsAsOne',true,'ReturnOnError',false);
    NoHeader{i} = [dataArray{1,i}{1,1} dataArray{1,i}{1,2} dataArray{1,i}{1,3}];
    fclose(fileID{i});
    
    clearvars fileID dataArray
    
end

% Reset origin

% % Initially set to these values to find the offsets used to center the data.
% xorigin = min(NoHeader{i}(:,1)); 
% yorigin = min(NoHeader{i}(:,2));
% xlength = max(NoHeader{i}(:,1))-min(NoHeader{i}(:,1));
% ylength = max(NoHeader{i}(:,2))-min(NoHeader{i}(:,2));
% 
% for i = 1:EndRunFiles-StartRunFiles+1
%     Displaced{i}(:, 1) = NoHeader{i}(:, 1) - xorigin;
%     Displaced{i}(:, 2) = NoHeader{i}(:, 2) - yorigin;
%     Displaced{i}(:, 3) = NoHeader{i}(:, 3);
% end

xorigin = 249 + 22.5;
yorigin = 198 + 12;
xlength = 178.3 - 20.5;
ylength = 92.45 - 13;

for i = 1:EndRunFiles-StartRunFiles+1
    Displaced1{i}(:, 1) = NoHeader{i}(:, 1) - xorigin;
    Displaced1{i}(:, 2) = NoHeader{i}(:, 2) - yorigin;
    Displaced1{i}(:, 3) = NoHeader{i}(:, 3);
end

%%

% This is used if there is need for a rigid rotation around the z-axis.
% Sometimes there is a 1-2 degree misalignment of the sample on the CMM;
% this corrects that issue.

% http://www.mathworks.com/matlabcentral/fileexchange/30864-3d-rotation-about-shifted-axis/content/AxelRot.m

anglex1 = 2.6;
angley1 = 1.75;
anglex2 = 155.6;
angley2 = 3.5;

for i = 1:EndRunFiles-StartRunFiles+1;
    [Displaced2{i} tempR{i} tempt{i}] = AxelRot(transpose(Displaced1{i}),atan(-(angley2-angley1)/(anglex2-anglex1))*180/pi,transpose([0,0,1]),transpose([0,0,0]));
    Displaced3{i} = transpose(Displaced2{i});
end

%%

% Fits a plane through the original data set, will be used later to
% reorient the data set so they are all nominally flat, data inspected is
% slightly smaller to account for the magnet placed on the surface.

for i = 1:EndRunFiles-StartRunFiles+1;
    OriginalXYZ{i} = Displaced3{i};
    TF1{i} = OriginalXYZ{i}(:,1)<0 | OriginalXYZ{i}(:,1)>xlength | OriginalXYZ{i}(:,2)<0 | OriginalXYZ{i}(:,2)>ylength;
    OriginalXYZ{i}(TF1{i}, :) = [];
end

for i = 1:EndRunFiles-StartRunFiles+1;
    [normalfit{i},vectorfit{i},pointfit{i}] = affine_fit(OriginalXYZ{i});
end

%%

desirednormal = [0 0 1];

% 3D rigid body rotation
% http://www.engr.uvic.ca/~mech410/lectures/4_2_RotateArbi.pdf
for i = 1:EndRunFiles-StartRunFiles+1
    thetaR{i} = acos(dot(normalfit{i},desirednormal)/((((normalfit{i}(1))^2+(normalfit{i}(2))^2+(normalfit{i}(3))^2)^0.5)*(((desirednormal(1))^2+(desirednormal(2))^2+(desirednormal(3))^2)^0.5)));
    CrossLine{i} = cross(normalfit{i},desirednormal);
end

for i = 1:EndRunFiles-StartRunFiles+1
    D{i} = [1 0 0 -pointfit{i}(1); 0 1 0 -pointfit{i}(2); 0 0 1 -pointfit{i}(3); 0 0 0 1];
    L{i} = ((CrossLine{i}(1))^2+(CrossLine{i}(2))^2+(CrossLine{i}(3))^2)^(0.5);
    V{i} = ((CrossLine{i}(2))^2+(CrossLine{i}(3))^2)^(0.5);
    Rx{i} = [1 0 0 0; 0 CrossLine{i}(3)/V{i} -CrossLine{i}(2)/V{i} 0; 0 CrossLine{i}(2)/V{i} CrossLine{i}(3)/V{i} 0; 0 0 0 1];
    Ry{i} = [V{i}/L{i} 0 -CrossLine{i}(1)/L{i} 0; 0 1 0 0; CrossLine{i}(1)/L{i} 0 V{i}/L{i} 0; 0 0 0 1];
    Rz{i} = [cos(thetaR{i}) -sin(thetaR{i}) 0 0; sin(thetaR{i}) cos(thetaR{i}) 0 0; 0 0 1 0; 0 0 0 1];
end

%%
for i = 1:EndRunFiles-StartRunFiles+1
    AndrewOriginalXYZ{i} = vertcat(Displaced3{i}',ones(1,length(Displaced3{i}')));
    RotatedAndrewOriginalXYZ{i} = (inv(D{i})*inv(Rx{i})*inv(Ry{i})*Rz{i}*Ry{i}*Rx{i}*D{i}*AndrewOriginalXYZ{i})';
    pointCloudAdjustedfit{i} = RotatedAndrewOriginalXYZ{i}(:,1:3);
    [normalfit2{i},vectorfit2{i},pointfit2{i}] = affine_fit(RotatedAndrewOriginalXYZ{i}(:,1:3));
end

%%

% Rotation of point cloud to match the local x and y with the global x and
% y after the normal rotation has been performed.

for i = 1:EndRunFiles-StartRunFiles+1
    pointCloudAdjustedfitUnrotated{i} = pointCloudAdjustedfit{i};
    [minxvalue{i},minxind{i}]=min(pointCloudAdjustedfitUnrotated{i}(:,1));
    [minyvalue{i},minyind{i}]=min(pointCloudAdjustedfitUnrotated{i}(:,2));
    [maxxvalue{i},maxxind{i}]=max(pointCloudAdjustedfitUnrotated{i}(:,1));
    [maxyvalue{i},maxyind{i}]=max(pointCloudAdjustedfitUnrotated{i}(:,2));
    th{i} = -atan((pointCloudAdjustedfitUnrotated{i}(maxxind{i},2)-pointCloudAdjustedfitUnrotated{i}(minyind{i},2))/(pointCloudAdjustedfitUnrotated{i}(maxxind{i},1)-pointCloudAdjustedfitUnrotated{i}(minyind{i},1)));
    pointCloudAdjustedfitOrigined{i}(:,1) = pointCloudAdjustedfitUnrotated{i}(:,1)-minxvalue{i};
    pointCloudAdjustedfitOrigined{i}(:,2) = pointCloudAdjustedfitUnrotated{i}(:,2)-minyvalue{i};
    pointCloudAdjustedfitOrigined{i}(:,3) = pointCloudAdjustedfitUnrotated{i}(:,3);
    pointCloudAdjustedfit{i} = pointCloudAdjustedfitUnrotated{i}*[cos(th{i}) sin(th{i}) 0; -sin(th{i}) cos(th{i}) 0; 0 0 1];
    Displaced{i} = pointCloudAdjustedfit{i};
end

%%

% Values to trunctate the y range and values to create the mesh in x and y.

xShiftFactor = 0;
yShiftFactor = 0;

XModifierLower = 0 + xShiftFactor;
XModifierUpper = xlength + xShiftFactor;

YModifierLower = 0 + yShiftFactor;
YModifierUpper = ylength + yShiftFactor;

for i = 1:EndRunFiles-StartRunFiles+1
    tablingvaluex{i} = 0.05;
    tablingvaluey{i} = tablingvaluex{i};
end

%%

% % Run using these values before the origin is known.
% for i = 1:EndRunFiles-StartRunFiles+1
%     MinimumBx{i} = min(Displaced{i}(:, 1));
%     MaximumBx{i} = max(Displaced{i}(:, 1));
%     MinimumBy{i} = min(Displaced{i}(:, 2));
%     MaximumBy{i} = max(Displaced{i}(:, 2));
% end

% Slightly truncates the data set due to some elements dropping along the
% fringes. Run this after the origin is known.
for i = 1:EndRunFiles-StartRunFiles+1
    MinimumBx{i} = XModifierLower; %min(Displaced{i}(:, 1));
    MaximumBx{i} = XModifierUpper; %max(Displaced{i}(:, 1));
    MinimumBy{i} = YModifierLower; %min(Displaced{i}(:, 2));
    MaximumBy{i} = YModifierUpper; %max(Displaced{i}(:, 2));
end

% Creates a mesh running from the minimum x to the maximum x and the
% adjusted minimum y to the adjusted maximum y in steps of tablingvaluex
% and tablingvaluey, respectively.
for i = 1:EndRunFiles-StartRunFiles+1
    [Xg{i},Yg{i}] = meshgrid(MinimumBx{i}:tablingvaluex{i}:MaximumBx{i}, MinimumBy{i}:tablingvaluey{i}:MaximumBy{i});
end

% Interpolation function for the z-dimension, returns 'none' when the
% extrapolation function is used.
for i = 1:EndRunFiles-StartRunFiles+1
    Zqf{i} = scatteredInterpolant(Displaced{i}(:, 1),Displaced{i}(:, 2),Displaced{i}(:, 3),'linear','nearest');
end

% Creates an interpolated surface from the Xg and Yg mesh.
for i = 1:EndRunFiles-StartRunFiles+1
    sprintf(strcat('\n','InterpolatedSurface ',num2str(i)))
    Zq{i}   = Zqf{i}(Xg{i},Yg{i});
end

%% Clears some data from RAM

clear Displaced NoHeader

%% Path length measurements

zMinForPlots = min(Zq{i}(:));

% Create a new array corresponding with lines constant Y

for i = 1:EndRunFiles-StartRunFiles+1
    for j = 1:floor((MaximumBy{i}-MinimumBy{i})/(tablingvaluey{i}))%1:((floor(MaximumBy{i})-ceil(MinimumBy{i}))/(tablingvaluey{i}))
        XZAlongConstY{i,j} = [transpose(Xg{i}(j,1:end)) transpose(Yg{i}(j,1:end)) transpose(Zq{i}(j,1:end)-zMinForPlots)];
    end
end

% Input this new array into the euclideanDistance calculator to measure the
% path length along lines of constant Y

for i = 1:EndRunFiles-StartRunFiles+1
    XZAlongConstYDistanceTotal{i} = zeros(length(1:((floor(MaximumBy{i})-ceil(MinimumBy{i}))/(tablingvaluey{i}))),3);
    for j = 1:floor((MaximumBy{i}-MinimumBy{i})/(tablingvaluey{i}))%1:((floor(MaximumBy{i})-ceil(MinimumBy{i}))/(tablingvaluey{i}))
        XZAlongConstYDistance{i,j} = [XZAlongConstY{i,j} EuclideanDistance(XZAlongConstY{i,j}(:, [1,3]))];
        XZAlongConstYDistanceTotal{i}(j,:) = [XZAlongConstYDistance{i,j}(1,2) XZAlongConstYDistance{i,j}(end,end) size(XZAlongConstYDistance{i,j},1)]; % ['y position' 'path length' 'number of elements']
    end
end

% Derivatives of z with respect to x to find the points needed to remove
% the table then polycarbonate backing

for i = 1:EndRunFiles-StartRunFiles+1
    for j = 1:floor((MaximumBy{i}-MinimumBy{i})/(tablingvaluey{i}))
        dZqbydXg{i,j} = [transpose(Xg{i}(j,:)) [0;transpose(rdivide(Zq{i}(j,2:end)-Zq{i}(j,1:end-1),Xg{i}(j,2:end)-Xg{i}(j,1:end-1)))]];
    end
end

for i = 1:EndRunFiles-StartRunFiles+1
    for j = 1:floor((MaximumBy{i}-MinimumBy{i})/(tablingvaluey{i}))
        d2ZqbydX2g{i,j} = [transpose(Xg{i}(j,:)) [0;transpose(rdivide(transpose(dZqbydXg{i,j}(2:end,2)-dZqbydXg{i,j}(1:end-1,2)),Xg{i}(j,2:end)-Xg{i}(j,1:end-1)))]];
    end
end

% Derivatives of z with respect to y to find regions away from edges to
% remove edge effects.
for i = 1:EndRunFiles-StartRunFiles+1
    for k = 1:floor((MaximumBx{i}-MinimumBx{i})/(tablingvaluex{i}))
        dZqbydYg{i,k} = [Yg{i}(:,1) [0;rdivide(Zq{i}(2:end,k)-Zq{i}(1:end-1,k),Yg{i}(2:end,k)-Yg{i}(1:end-1,k))]];
    end
end

%%

% Finds the first y position such that the edges have been removed (ie.
% where the surface is just prepreg and not polycarbonate).

dataclippingminy = 0;
dataclippingmaxy = ylength;

polycarbonatePrepregLeft = 10;
polycarbonatePrepregRight = 150;

regionsY = [0 10; 74 85; 150 xlength];

minPeakDistanceYZConstX = 1;

filterFrameLength = 21;

for i = 1:EndRunFiles-StartRunFiles+1;

    FirstYValFinderArray{i} = [];
    LastYValFinderArray{i} = [];
    
    FirstYValFinderArray1{i} = [];
    LastYValFinderArray1{i} = [];
    
    for region = 1:length(regionsY)-1
    
        for k = ceil(((regionsY(region,2)-MinimumBx{i})/(tablingvaluex{i}))):floor(((regionsY(region+1,1)-MinimumBx{i})/(tablingvaluex{i})));

            smoothedYZConstX{i,k} = sgolayfilt(horzcat(Yg{i}(:,k),Zq{i}(:,k)),3,41);
            dzdyDataSmoothedYZConstX{i,k} = [smoothedYZConstX{i,k}(:,1) [0;diff(smoothedYZConstX{i,k}(:,2))./diff(smoothedYZConstX{i,k}(:,1))]];

            lowerClippingMaskYZConstX{i,k} = (smoothedYZConstX{i,k}(:,1) > dataclippingminy) & (smoothedYZConstX{i,k}(:,1) < (max(Yg{i}(1:end,k))-min(Yg{i}(1:end,k)))*0.2);
            higherClippingMaskYZConstX{i,k} = (smoothedYZConstX{i,k}(:,1) > (max(Yg{i}(1:end,k))-min(Yg{i}(1:end,k)))*0.8) & (smoothedYZConstX{i,k}(:,1) < dataclippingmaxy);

            smoothedLowerYZConstX{i,k} = sgolayfilt(smoothedYZConstX{i,k}(lowerClippingMaskYZConstX{i,k}, :),3,filterFrameLength);
            smoothedHigherYZConstX{i,k} = sgolayfilt(smoothedYZConstX{i,k}(higherClippingMaskYZConstX{i,k}, :),3,filterFrameLength);
            smoothedLowerdYdZConstX{i,k} = sgolayfilt(dzdyDataSmoothedYZConstX{i,k}(lowerClippingMaskYZConstX{i,k}, :),3,filterFrameLength);
            smoothedHigherdYdZConstX{i,k} = sgolayfilt(dzdyDataSmoothedYZConstX{i,k}(higherClippingMaskYZConstX{i,k}, :),3,filterFrameLength);

            [pospksSmoothedLowerYZConstX{i,k} poslocsSmoothedLowerYZConstX{i,k}] = findpeaks(smoothedLowerYZConstX{i,k}(:,2),smoothedLowerYZConstX{i,k}(:,1),'MinPeakDistance',minPeakDistanceYZConstX);

            [pospksSmoothedLowerdYdZConstX{i,k} poslocsSmoothedLowerdYdZConstX{i,k}] = findpeaks(smoothedLowerdYdZConstX{i,k}(:,2),smoothedLowerdYdZConstX{i,k}(:,1),'MinPeakDistance',minPeakDistanceYZConstX);
            [negpksSmoothedLowerdYdZConstX{i,k} neglocsSmoothedLowerdYdZConstX{i,k}] = findpeaks(-smoothedLowerdYdZConstX{i,k}(:,2),smoothedLowerdYdZConstX{i,k}(:,1),'MinPeakDistance',minPeakDistanceYZConstX);

            [pospksSmoothedHigherYZConstX{i,k} poslocsSmoothedHigherYZConstX{i,k}] = findpeaks(smoothedHigherYZConstX{i,k}(:,2),smoothedHigherYZConstX{i,k}(:,1),'MinPeakDistance',minPeakDistanceYZConstX);

            [pospksSmoothedHigherdYdZConstX{i,k} poslocsSmoothedHigherdYdZConstX{i,k}] = findpeaks(smoothedHigherdYdZConstX{i,k}(:,2),smoothedHigherdYdZConstX{i,k}(:,1),'MinPeakDistance',minPeakDistanceYZConstX);
            [negpksSmoothedHigherdYdZConstX{i,k} neglocsSmoothedHigherdYdZConstX{i,k}] = findpeaks(-smoothedHigherdYdZConstX{i,k}(:,2),smoothedHigherdYdZConstX{i,k}(:,1),'MinPeakDistance',minPeakDistanceYZConstX);

            [valMaxPeakSmoothedLowerdYdZConstX{i,k} posMaxPeakSmoothedLowerdYdZConstX{i,k}] = max(pospksSmoothedLowerdYdZConstX{i,k});
            [valMaxPeakSmoothedHigherdYdZConstX{i,k} posMaxPeakSmoothedHigherdYdZConstX{i,k}] = max(negpksSmoothedHigherdYdZConstX{i,k});

            greaterPeaksSmoothLowerYZConstX{i,k} = min(poslocsSmoothedLowerYZConstX{i,k}(poslocsSmoothedLowerYZConstX{i,k}-poslocsSmoothedLowerdYdZConstX{i,k}(pospksSmoothedLowerdYdZConstX{i,k}==max(pospksSmoothedLowerdYdZConstX{i,k}))>0));
            greaterPeaksSmoothHigherYZConstX{i,j} = max(poslocsSmoothedHigherYZConstX{i,k}((poslocsSmoothedHigherYZConstX{i,k}-neglocsSmoothedHigherdYdZConstX{i,k}(negpksSmoothedHigherdYdZConstX{i,k}==-min(-negpksSmoothedHigherdYdZConstX{i,k}))<0)));

            if length(greaterPeaksSmoothLowerYZConstX{i,k}) >=1
                FirstYValFinder{i,k} = greaterPeaksSmoothLowerYZConstX{i,k}(1);
            else
                FirstYValFinder{i,k} = NaN;
            end

            if length(greaterPeaksSmoothHigherYZConstX{i,j}) >=1
                LastYValFinder{i,k} = greaterPeaksSmoothHigherYZConstX{i,j};
            else
                LastYValFinder{i,k} = NaN;
            end

            FirstYValFinderArray{i} = [FirstYValFinderArray{i};k*tablingvaluex{i}+MinimumBx{i} FirstYValFinder{i,k}];
            LastYValFinderArray{i} = [LastYValFinderArray{i};k*tablingvaluex{i}+MinimumBx{i} LastYValFinder{i,k}];

        end
        
    end
end

%%

clearvars coeffsforFirstYVal coeffsforLastYVal FirstYVal LastYVal

for i = 1:EndRunFiles-StartRunFiles+1;
    coeffsforFirstYVal{i} = polyfit(FirstYValFinderArray{i}(:,1), FirstYValFinderArray{i}(:,2), 1);
    coeffsforLastYVal{i} = polyfit(LastYValFinderArray{i}(:,1), LastYValFinderArray{i}(:,2), 1);
    
    FirstYVal{i} = max(polyval(coeffsforFirstYVal{i},MinimumBx{i}),polyval(coeffsforFirstYVal{i},MaximumBx{i}));
    LastYVal{i} = min(polyval(coeffsforLastYVal{i},MinimumBx{i}),polyval(coeffsforLastYVal{i},MaximumBx{i}));
end

%%

% Standard deviation and mean of the first derivative after x clipping

dataclippingminx1 = 0;
dataclippingmaxx1 = 6;
dataclippingminx2 = 151;
dataclippingmaxx2 = xlength;

minPeakDistanceXZConstY = 0.75;

filterFrameLength = 21; % Higher filtering value increases the 'blurring', ie. higher effectiveness (less errors due to short disturbances in the first derivative) but decreases the resolution.

for i = 1:EndRunFiles-StartRunFiles+1
    
    FirstXValFinderArray1{i} = [];
    LastXValFinderArray1{i} = [];
    
    firstXArray{i} = [];
    lastXArray{i} = [];
    
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        
        smoothedXZConstY{i,j} = sgolayfilt(horzcat(Xg{i}(j,:)',Zq{i}(j,:)'),3,filterFrameLength);
        dzdxDataSmoothedXZConstY{i,j} = [smoothedXZConstY{i,j}(:,1) [0;diff(smoothedXZConstY{i,j}(:,2))./diff(smoothedXZConstY{i,j}(:,1))]];
        
        lowerClippingMaskXZConstY{i,j} = (smoothedXZConstY{i,j}(:,1) > dataclippingminx1) & (smoothedXZConstY{i,j}(:,1) < dataclippingmaxx1);
        higherClippingMaskXZConstY{i,j} = (smoothedXZConstY{i,j}(:,1) > dataclippingminx2) & (smoothedXZConstY{i,j}(:,1) < dataclippingmaxx2);

        smoothedLowerXZConstY{i,j} = smoothedXZConstY{i,j}(lowerClippingMaskXZConstY{i,j}, :);
        smoothedHigherXZConstY{i,j} = smoothedXZConstY{i,j}(higherClippingMaskXZConstY{i,j}, :);
        
        smoothedLowerdXdZConstY{i,j} = sgolayfilt(dzdxDataSmoothedXZConstY{i,j}(lowerClippingMaskXZConstY{i,j}, :),3,filterFrameLength);
        smoothedHigherdXdZConstY{i,j} = sgolayfilt(dzdxDataSmoothedXZConstY{i,j}(higherClippingMaskXZConstY{i,j}, :),3,filterFrameLength);

        FirstX{i,j} = 3.75; % Manually set based on the figure outputs
        LastX{i,j} = 155; % Manually set based on the figure outputs
        
        firstXArray{i} = [firstXArray{i}; FirstX{i,j} j*tablingvaluey{i}+MinimumBy{i} Zq{i}(j,floor((FirstX{i,j}+MinimumBx{i})/tablingvaluex{i}))]; % ['x' 'y' 'z'] of first x position at y location
        lastXArray{i} = [lastXArray{i}; LastX{i,j} j*tablingvaluey{i}+MinimumBy{i} Zq{i}(j,floor((LastX{i,j}+MinimumBx{i})/tablingvaluex{i}))];
        
        clearvars lowerClippingMaskXZConstY higherClippingMaskXZConstY smoothedLowerdXdZConstY...
            smoothedHigherdXdZConstY pospksSmoothedLowerdXdZConstY poslocsSmoothedLowerdXdZConstY negpksSmoothedLowerdXdZConstY neglocsSmoothedLowerdXdZConstY...
            pospksSmoothedHigherdXdZConstY poslocsSmoothedHigherdXdZConstY negpksSmoothedHigherdXdZConstY neglocsSmoothedHigherdXdZConstY...
            valMaxPeakSmoothedLowerdXdZConstY posMaxPeakSmoothedLowerdXdZConstY valMaxPeakSmoothedHigherdXdZConstY posMaxPeakSmoothedHigherdXdZConstY...
            greaterPeaksSmootherLowerdXdZConstY lesserPeaksSmootherHigherdXdZConstYTEMP1 lesserPeaksSmootherHigherdXdZConstYTEMP2 lesserPeaksSmootherHigherdXdZConstY
        
    end
    
end

%%

% Plots the smoothed beginning of the line of constant y, the first
% derivative, then the position corresponding with the first x value
% chosen.  This is done by finding the minimum peak to the left or right,
% depending on if the maximum or minimum value is of interest, of the
% maximum peak of the first derivative.

% i = 1;
% 
% for j = 267;%ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})-ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}))/10):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
%     figure;
%         plotyy(smoothedXZConstY{i,j}(lowerClippingMaskXZConstY{i,j},1),smoothedXZConstY{i,j}(lowerClippingMaskXZConstY{i,j},2),smoothedLowerdXdZConstY{i,j}(:,1),smoothedLowerdXdZConstY{i,j}(:,2));
%         hold on
%         y1=get(gca,'ylim');
%         hold on
%         plot([FirstXValFinder{i,j} FirstXValFinder{i,j}],y1);
% end

%%

% Similar to above, except corresponding to the right most x value.

% i=1;
% 
% for j = 266;%ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})-ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}))/10):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
%     figure;
%         yyaxis left
%         plot(smoothedXZConstY{i,j}(higherClippingMaskXZConstY{i,j},1),smoothedXZConstY{i,j}(higherClippingMaskXZConstY{i,j},2));
%         hold on
%         
%         scatter(maxpeaksArrayHigherXZConstY{i,j}(:,1),maxpeaksArrayHigherXZConstY{i,j}(:,2),'wv','MarkerFaceColor','r')
%         scatter(minpeaksArrayHigherXZConstY{i,j}(:,1),minpeaksArrayHigherXZConstY{i,j}(:,2),'wv','MarkerFaceColor','b')
%         hold off
%         
%         yyaxis right
%         plot(smoothedHigherdXdZConstY{i,j}(:,1),smoothedHigherdXdZConstY{i,j}(:,2));
%         hold on
%         scatter(neglocsSmoothedHigherdXdZConstY{i,j},-negpksSmoothedHigherdXdZConstY{i,j})
%         
%         title(num2str(j))
%         
%         hold off
%         y1=get(gca,'ylim');
%         hold on
%         plot([LastXValFinder{i,j} LastXValFinder{i,j}],y1);
% end

%%

clearvars XgvsZq inRangeXgvsZqLower inRangeXgvsZqUpper inRangeXgvsZq circleParameters...
    circleParametersIndividual signCircleIndividual meanCurvature stdCurvature

% Determines the x-z curvature for each frame and line of constant y

% Creates the range of values to investigate, removes the wrinkle using the
% clipping used earlier.

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/(tablingvaluey{i})):floor((LastYVal{i}-MinimumBy{i})/(tablingvaluey{i}))
        XgvsZq{i,j} = [transpose(Xg{i}(j,:)) transpose(Zq{i}(j,:))];
        
        inRangeXgvsZq{i,j} = [];
        inRangeXgvsZqRegion{i,j,1} = XgvsZq{i,j}(XgvsZq{i,j}(:,1)>=FirstX{i,j} & XgvsZq{i,j}(:,1)<=regionsY(1,2),:);
        inRangeXgvsZqRegion{i,j,length(regionsY)} = XgvsZq{i,j}(XgvsZq{i,j}(:,1)>=regionsY(length(regionsY),1) & XgvsZq{i,j}(:,1)<=LastX{i,j},:);
        inRangeXgvsZq{i,j} = [inRangeXgvsZq{i,j};inRangeXgvsZqRegion{i,j,1}];
        for region = 2:length(regionsY)-1
            inRangeXgvsZqRegion{i,j,region} = XgvsZq{i,j}(XgvsZq{i,j}(:,1)>=regionsY(region-1,2) & XgvsZq{i,j}(:,1)<=regionsY(region,1),:);
            
            inRangeXgvsZq{i,j} = [inRangeXgvsZq{i,j};inRangeXgvsZqRegion{i,j,region}];
        end
        inRangeXgvsZq{i,j} = [inRangeXgvsZq{i,j};inRangeXgvsZqRegion{i,j,length(regionsY)}];
        
    end
end

for i = 1:EndRunFiles-StartRunFiles+1
    circleParameters{i} = [];
    for j = ceil((FirstYVal{i}-MinimumBy{i})/(tablingvaluey{i})):floor((LastYVal{i}-MinimumBy{i})/(tablingvaluey{i}))
        circleParametersIndividual{i,j} = CircleFit(inRangeXgvsZq{i,j});
        signCircleIndividual{i,j} = circleParametersIndividual{i,j}(2)>Zqf{i}(circleParametersIndividual{i,j}(1),(j-1)*tablingvaluey{i}+MinimumBy{i}); % If the Z value of the real curve is less than the Z value of the fitted circle at the same centroid x value, then it is very likely that the curvature is negative, meaning a smiley face shape.
        circleParameters{i} = [circleParameters{i};(j-1)*tablingvaluey{i}+MinimumBy{i} 1/circleParametersIndividual{i,j}(3) signCircleIndividual{i,j}];
    end
    meanCurvature{i} = mean(circleParameters{i}(:,2));
    stdCurvature{i} = nanstd(circleParameters{i}(:,2));
%     signCurvature{i} = sum(circleParameters{i}(:,3)<0.5); %  If this value is non zero, then there are some lines of constant y that have positive curvatures.
end

%%

% Plots the fitted circle overlayed on the raw data.  This is to visualize
% the curvature direction as well as to see if the fit is reasonable.

% i = 1;
% j = 800;
% 
% figure;
% plot(Xg{i}(j,:),Zq{i}(j,:));
% hold on
% CirclePlot(circleParametersIndividual{i,j}(1),circleParametersIndividual{i,j}(2),circleParametersIndividual{i,j}(3),10000); %10000 corresponds to a fineness of the circle output, lower = lower resolution circle on plot
% axis([MinimumBx{i} MaximumBx{i} min(inRangeXgvsZq{i,j}(:,2))-0.2 max(inRangeXgvsZq{i,j}(:,2))+0.2])

%%

clear FirstXValFinderArray1 LastXValFinderArray1 dzdxDataSmoothedXZConstY lowerClippingMaskXZConstY higherClippingMaskXZConstY...
    smoothedLowerdXdZConstY smoothedHigherdXdZConstY pospksSmoothedLowerdXdZConstY poslocsSmoothedLowerdXdZConstY negpksSmoothedLowerdXdZConstY...
    neglocsSmoothedLowerdXdZConstY pospksSmoothedHigherdXdZConstY poslocsSmoothedHigherdXdZConstY negpksSmoothedHigherdXdZConstY neglocsSmoothedHigherdXdZConstY...
    valMaxPeakSmoothedLowerdXdZConstY posMaxPeakSmoothedLowerdXdZConstY valMaxPeakSmoothedHigherdXdZConstY posMaxPeakSmoothedHigherdXdZConstY{i,j}...
    greaterPeaksSmootherLowerdXdZConstY lesserPeaksSmootherHigherdXdZConstYTEMP1 ...
    lesserPeaksSmootherHigherdXdZConstYTEMP2 lesserPeaksSmootherHigherdXdZConstY FirstXValFinder...
    LastXValFinder FirstXValFinderArray1 LastXValFinderArray1

%%
% Creates the appropriate clipping mask for each increment in regular xz
% space

% for i = 1:EndRunFiles-StartRunFiles+1
%     for j = 1:floor((MaximumBy{i}-MinimumBy{i})/(tablingvaluey{i}))
%         RealSpaceClipping{i,j} = transpose(Xg{i}(j,1:end)>FirstX{i,j}&Xg{i}(j,1:end)<LastX{i,j});
%         PreClippedXg{i,j} = Xg{i}(j,1:end);
%         ClippedXg{i,j} = PreClippedXg{i}(RealSpaceClipping{i,j});
%         PreClippedZq{i,j} = Zq{i}(j,1:end);
%         ClippedZq{i,j} = PreClippedZq{i,j}(RealSpaceClipping{i,j});
%     end
% end

%%

% clear RealSpaceClipping PreClippedXg PreClippedZq

%%

% % Measures the lengths of the clipped lines
% 
% for i = 1:EndRunFiles-StartRunFiles+1
%     ClippedXZAlongConstYDistanceTotal{i} = zeros(length(1:((floor(MaximumBy{i})-ceil(MinimumBy{i}))/(tablingvaluey{i}))),3);
%     for j = 1:floor((MaximumBy{i}-MinimumBy{i})/(tablingvaluey{i}))%1:((floor(MaximumBy{i})-ceil(MinimumBy{i}))/(tablingvaluey{i}))
%         ClippedXZAlongConstYDistance{i,j} = [transpose(ClippedXg{i,j}) transpose(repmat((j-1)*tablingvaluey{i}+MinimumBy{i},1,length(ClippedXg{i,j}))) transpose(ClippedZq{i,j}) EuclideanDistance([transpose(ClippedXg{i,j}) transpose(ClippedZq{i,j})])];
%         ClippedXZAlongConstYDistanceTotal{i}(j,:) = [ClippedXZAlongConstYDistance{i,j}(1,2) ClippedXZAlongConstYDistance{i,j}(end,end) size(ClippedXZAlongConstYDistance{i,j},1)]; % ['y position' 'path length' 'number of elements']
%     end
% end

%%

% clear ClippedXg ClippedZq

% Outcome, the path length changes along Y with what appears to be a
% similar change in the number of elements being counted.  Solution,
% resample the clipped space so that the number of elements being counted
% along each Y is the same

%%

% Take an average for all frames, the point here needs to be constant
ClippedElementCount = 6500;%(200-0)/(tablingvaluey{1});

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        ConstElementClippedXg{i,j} = FirstX{i,j}:(LastX{i,j}-FirstX{i,j})/(ClippedElementCount-1):LastX{i,j};
        ConstElementClippedZq{i,j} = Zqf{i}(ConstElementClippedXg{i,j},repmat((j-1)*tablingvaluey{i}+MinimumBy{i},1,length(ConstElementClippedXg{i,j})));
        ConstElementClippedXgZq{i,j} = [transpose(ConstElementClippedXg{i,j}) transpose(ConstElementClippedZq{i,j})];
        
        clearvars ConstElementClippedZq
    end
end

%%

clearvars Zqf

%%

% Measures the lengths of the clipped lines with now constant number of
% elements

for i = 1:EndRunFiles-StartRunFiles+1
    ConstElementClippedXZAlongConstYDistanceTotal{i} = zeros(length(1:((floor(MaximumBy{i})-ceil(MinimumBy{i}))/(tablingvaluey{i}))),3);
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        ConstElementClippedXZAlongConstYDistance{i,j} = [ConstElementClippedXgZq{i,j}(:,1) repmat((j-1)*tablingvaluey{i}+MinimumBy{i},length(ConstElementClippedXg{i,j}),1) ConstElementClippedXgZq{i,j}(:,2) EuclideanDistance(ConstElementClippedXgZq{i,j})];
        ConstElementClippedXZAlongConstYDistanceTotal{i}(j,:) = [ConstElementClippedXZAlongConstYDistance{i,j}(1,2) ConstElementClippedXZAlongConstYDistance{i,j}(end,end) size(ConstElementClippedXZAlongConstYDistance{i,j},1)]; % ['y position' 'path length' 'number of elements']
    end
    
end

clearvars ConstElementClippedXg ConstElementClippedXgZq

%%

clearvars pzx szx muzx f_zx DetrendedXZAlongConstYDistance

% 2016-04-16 Detrending may be required to find a more consistent method of
% determining the minimum peak values.  Ignore the next statement for the
% time being.
% Detrending didn't work well.

% Detrending the data for peak detection.
% http://www.mathworks.com/help/signal/examples/peak-analysis.html

polynomialfitdegree = 4;

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        [pzx{i,j},szx{i,j},muzx{i,j}] = polyfit(ConstElementClippedXZAlongConstYDistance{i,j}(:,1),ConstElementClippedXZAlongConstYDistance{i,j}(:,3),polynomialfitdegree);
        % Least-square fits the data to a polynomial of degree
        % polynomialfitdegree
        % http://www.mathworks.com/help/matlab/ref/polyfit.html
        f_zx{i,j} = polyval(pzx{i,j},ConstElementClippedXZAlongConstYDistance{i,j}(:,1),[],muzx{i,j});
        % Evaluates the fitted function at each point of the original x set
        DetrendedXZAlongConstYDistance{i,j} = ConstElementClippedXZAlongConstYDistance{i,j};
        DetrendedXZAlongConstYDistance{i,j}(:,3) = ConstElementClippedXZAlongConstYDistance{i,j}(:,3)-f_zx{i,j};
        % Removes the trend data from the original data set
        
        clearvars szx pzx muzx f_zx
        
    end
end

%%

% i=1;
% j=500;
% 
% figure;
% plot(DetrendedXZAlongConstYDistance{i,j}(:,1),DetrendedXZAlongConstYDistance{i,j}(:,3),ConstElementClippedXZAlongConstYDistance{i,j}(:,1),f_zx{i,j}-mean(ConstElementClippedXZAlongConstYDistance{i,j}(:,3)),XZAlongConstYDistance{i,j}(:,1),XZAlongConstYDistance{i,j}(:,3)-mean(XZAlongConstYDistance{i,j}(:,3)))

%%

clearvars SmoothedConstElementClippedXZAlongConstYDistance

% Smooth the data to remove arbitrary peaks
% http://www.mathworks.com/help/signal/ref/sgolayfilt.html

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        SmoothedConstElementClippedXZAlongConstYDistance{i,j} = sgolayfilt(DetrendedXZAlongConstYDistance{i,j}(:,[1,3]),3,101);
        nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j} = sgolayfilt(ConstElementClippedXZAlongConstYDistance{i,j}(:,[1,3]),3,101);
    end
end

%%

clearvars pks_MaxPeakXZConstY locs_MaxPeakXZConstY pks_MinPeakXZConstYproximity locs_MinPeakXZConstYproximity minpeakdistanceXZConstYmax minpeakdistanceXZConstYminminpeakheightXZConstYmin prominenceMinPeakVal

% Peak detection
% http://www.mathworks.com/help/signal/examples/peak-analysis.html

minpeakdistanceXZConstYmax = 2;
minpeakdistanceXZConstYmin = 1;

prominenceMaxPeakVal = 0.001;
prominenceMinPeakVal = 0.001; % This value is a tuning value.  It will cause the most headaches. Tune it based on the outputs from the graph a few sections below. http://www.mathworks.com/help/signal/ref/findpeaks.html

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        [pks_MaxPeakXZConstY{i,j},locs_MaxPeakXZConstY{i,j}] = findpeaks(SmoothedConstElementClippedXZAlongConstYDistance{i,j}(:,2),SmoothedConstElementClippedXZAlongConstYDistance{i,j}(:,1),'MinPeakDistance',minpeakdistanceXZConstYmax,'MinPeakProminence',prominenceMaxPeakVal);
        [pks_MinPeakXZConstY{i,j},locs_MinPeakXZConstY{i,j}] = findpeaks(-SmoothedConstElementClippedXZAlongConstYDistance{i,j}(:,2),SmoothedConstElementClippedXZAlongConstYDistance{i,j}(:,1),'MinPeakDistance',minpeakdistanceXZConstYmin,'MinPeakProminence',prominenceMinPeakVal);
    end
end

%%

% Finds the z-peak corresponding with the x-location in the (original)
% non-detrended data set using the peak values determined in the
% detrended data set.

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        for k = 1:length(locs_MaxPeakXZConstY{i,j})
            peakZNonDetrended{i,j,k} = nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(:,1)==locs_MaxPeakXZConstY{i,j}(k),2);
        end
    end
end

%%

clearvars closestRightPeakXTemp closestLeftPeakXTemp closestRightPeakX closestLeftPeakX closestLeftPeakXInd closestRightPeakXInd closestLeftPeakZ closestRightPeakZ

%%

% Determines the two nearest minimum peaks to each maximum peak, both in X
% and in Z.  This will be used to calculate the wrinkle height, and width,
% of each wrinkle according to the local maximum peak.  This will then be
% used to determine the appropriate value for the maximum wrinkle.

minimumDistanceBetweenMaxAndMinPeak = 1;

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        for k = 1:length(locs_MaxPeakXZConstY{i,j});

            closestRightPeakXTemp{i,j,k} = locs_MinPeakXZConstY{i,j}(locs_MinPeakXZConstY{i,j}>locs_MaxPeakXZConstY{i,j}(k)+minimumDistanceBetweenMaxAndMinPeak);
            closestLeftPeakXTemp{i,j,k} = locs_MinPeakXZConstY{i,j}(locs_MinPeakXZConstY{i,j}<locs_MaxPeakXZConstY{i,j}(k)-minimumDistanceBetweenMaxAndMinPeak);
            
            if length(closestRightPeakXTemp{i,j,k})<1 && length(closestLeftPeakXTemp{i,j,k})<1
                closestRightPeakX{i,j,k} = locs_MaxPeakXZConstY{i,j}(k);
                closestLeftPeakX{i,j,k} = locs_MaxPeakXZConstY{i,j}(k);
                
                closestRightPeakZ{i,j,k} = locs_MaxPeakXZConstY{i,j}(k);
                closestLeftPeakZ{i,j,k} = locs_MaxPeakXZConstY{i,j}(k);
                
            elseif length(closestRightPeakXTemp{i,j,k})<1
                closestRightPeakX{i,j,k} = locs_MaxPeakXZConstY{i,j}(k);
                closestLeftPeakX{i,j,k} = closestLeftPeakXTemp{i,j,k}(end);
                
                closestLeftPeakXInd{i,j,k} = find(locs_MinPeakXZConstY{i,j}==closestLeftPeakX{i,j,k},1);
                
                closestRightPeakZ{i,j,k} = locs_MaxPeakXZConstY{i,j}(k);
                closestLeftPeakZ{i,j,k} = -pks_MinPeakXZConstY{i,j}(closestLeftPeakXInd{i,j,k});
                
            elseif length(closestLeftPeakXTemp{i,j,k})<1
                closestLeftPeakX{i,j,k} = locs_MaxPeakXZConstY{i,j}(k);
                closestRightPeakX{i,j,k} = closestRightPeakXTemp{i,j,k}(1);
                
                closestRightPeakXInd{i,j,k} = find(locs_MinPeakXZConstY{i,j}==closestRightPeakX{i,j,k},1);
                
                closestLeftPeakZ{i,j,k} = locs_MaxPeakXZConstY{i,j}(k);
                closestRightPeakZ{i,j,k} = -pks_MinPeakXZConstY{i,j}(closestRightPeakXInd{i,j,k});
                
            else
                closestRightPeakX{i,j,k} = closestRightPeakXTemp{i,j,k}(1);
                closestLeftPeakX{i,j,k} = closestLeftPeakXTemp{i,j,k}(end);
                
                closestLeftPeakXInd{i,j,k} = find(locs_MinPeakXZConstY{i,j}==closestLeftPeakX{i,j,k},1);
                closestRightPeakXInd{i,j,k} = find(locs_MinPeakXZConstY{i,j}==closestRightPeakX{i,j,k},1);
            
                % Need the following two peak variables for the plot,
                % not used for anything else.
                closestLeftPeakZDetrended{i,j,k} = -pks_MinPeakXZConstY{i,j}(closestLeftPeakXInd{i,j,k});
                closestRightPeakZDetrended{i,j,k} = -pks_MinPeakXZConstY{i,j}(closestRightPeakXInd{i,j,k});
                
                closestLeftPeakZ{i,j,k} = nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(:,1)==closestLeftPeakX{i,j,k},2);
                closestRightPeakZ{i,j,k} = nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(:,1)==closestRightPeakX{i,j,k},2);
                
            end

        end
    end
end

%%

clearvars closestLeftPeakXInd closestRightPeakXInd closestRightPeakXTemp closestLeftPeakXTemp

%%

clear wrinkleHeight wrinkleArrayWithSmall wrinkleArray sortedWrinkleArray

% Determines the height and width of each wrinkle.

% Then removes all wrinkles which aren't at least realProminenceMaxPeakVal
% high.  I think the findpeaks function doesn't necessarily choose the real
% data max prominence so it's easier to set a very low value above and then
% manually only choose the 'tall' wrinkles.

realProminenceMaxPeakVal = 0.0625;

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        for k = 1:length(locs_MaxPeakXZConstY{i,j});
            wrinklePathLength{i,j,k} = EuclideanDistanceReturnScalar(nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(:,1)>=closestLeftPeakX{i,j,k} & nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(:,1)<=closestRightPeakX{i,j,k},:));
            wrinkleHeight{i,j,k} = [peakZNonDetrended{i,j,k}-(closestLeftPeakZ{i,j,k}+closestRightPeakZ{i,j,k})/2 sqrt((closestLeftPeakX{i,j,k}-closestRightPeakX{i,j,k}).^2+(closestLeftPeakZ{i,j,k}-closestRightPeakZ{i,j,k}).^2) closestLeftPeakX{i,j,k} closestRightPeakX{i,j,k} locs_MaxPeakXZConstY{i,j}(k) j*tablingvaluey{i}+MinimumBy{i} peakZNonDetrended{i,j,k} wrinklePathLength{i,j,k}]; %[ wrinkle height, wrinkle width, left x position, right x position, max height x position, max height y position, max height z position, wrinkle excess length]
%             wrinkleWidth{i,j,k} = sqrt((closestLeftPeakX{i,j,k}-closestRightPeakX{i,j,k}).^2+(closestLeftPeakZ{i,j,k}-closestRightPeakZ{i,j,k}).^2);
        end
        wrinkleArrayWithSmall{i,j} = vertcat(wrinkleHeight{i,j,1:length(locs_MaxPeakXZConstY{i,j})});
        wrinkleArray{i,j} = wrinkleArrayWithSmall{i,j}(wrinkleArrayWithSmall{i,j}(:,1)>=realProminenceMaxPeakVal,:);
        sortedWrinkleArray{i,j} = sortrows(wrinkleArray{i,j});
    end
end

%%

clearvars nWrinkles wrinkleArrayCondensed

% Condenses the top n wrinkles from each line of constant y into a single
% array.

% nWrinkles = 3;

wrinkleArrayCondensed{i} = [];

for i = 1:EndRunFiles-StartRunFiles+1
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
%         wrinkleArrayCondensed{i} = [wrinkleArrayCondensed{i};sortedWrinkleArray{i,j}(end-(nWrinkles-1):end,5:7)];
        wrinkleArrayCondensed{i} = [wrinkleArrayCondensed{i};sortedWrinkleArray{i,j}(:,5:7) sortedWrinkleArray{i,j}(:,1:4) sortedWrinkleArray{i,j}(:,8)];
    end
end

%%

% This is here because some frames don't have any features larger than the
% minimum feature size, therefore any math applied on those frames will
% error.   This is to correct that.

clearvars framesWithWrinkles

framesWithWrinkles = [];

for i = 1:EndRunFiles-StartRunFiles+1
    if length(wrinkleArrayCondensed{i} > 0)
        framesWithWrinkles = [framesWithWrinkles i];
    end
end

%%

clearvars tempval indx indy wrinklesOnlyArray

% Creates a surface array with only the maximum wrinkles positions as real
% values.

for i = framesWithWrinkles
    wrinklesOnlyArray{i} = NaN(size(Zq{i}));
    for j = ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
        for k = 1:length(sortedWrinkleArray{i,j}(:,1))
            [tempval indx{i,j,k}] = min(abs(Xg{i}(1,:)-sortedWrinkleArray{i,j}(end-(k-1),5)));
            [tempval indy{i,j,k}] = min(abs(Yg{i}(:,1)-sortedWrinkleArray{i,j}(end-(k-1),6)));
            wrinklesOnlyArray{i}(indy{i,j,k},indx{i,j,k}) = Zq{i}(indy{i,j,k},indx{i,j,k});
        end
    end
end

%%

% Sets a grid size to overlay onto the wrinkle only surface array which
% returns true for the window of the wrinkle.  Effectively increases the
% wrinkled domain from a single point to a grid.  Multiple close by wrinkles
% will constructively interfere and allow for one wrinkle to be demarcated
% from multiple discrete wrinkle points.

clearvars xWindowSize yWindowSize

xWindowSize = 10;
yWindowSize = 15;

%%

clearvars logicalWindowedWrinklesOnlyArray indexValuesOfNonZero windowedWrinklesOnlyArray zzeroedWindowedWrinklesOnlyArray

for i = framesWithWrinkles
    logicalWindowedWrinklesOnlyArray{i} = false(size(meshgrid(MinimumBx{i}:tablingvaluex{i}:MaximumBx{i}, MinimumBy{i}:tablingvaluey{i}:MaximumBy{i})));
    windowedWrinklesOnlyArray{i} = NaN(size(meshgrid(MinimumBx{i}:tablingvaluex{i}:MaximumBx{i}, MinimumBy{i}:tablingvaluey{i}:MaximumBy{i})));
    zzeroedWindowedWrinklesOnlyArray{i} = NaN(size(meshgrid(MinimumBx{i}:tablingvaluex{i}:MaximumBx{i}, MinimumBy{i}:tablingvaluey{i}:MaximumBy{i})));
    for x = 1:floor(((MaximumBx{i}-MinimumBx{i})/tablingvaluex{i})/xWindowSize)
        for y = 1:floor(((MaximumBy{i}-MinimumBy{i})/tablingvaluey{i})/yWindowSize)
            logicalWindowedWrinklesOnlyArray{i}((y-1)*yWindowSize+1:y*yWindowSize,(x-1)*xWindowSize+1:x*xWindowSize) = nansum(nansum(wrinklesOnlyArray{i}((y-1)*yWindowSize+1:y*yWindowSize,(x-1)*xWindowSize+1:x*xWindowSize)));
%             indexValuesOfNonZero{i} = find(logicalWindowedWrinklesOnlyArray{i});
%             windowedWrinklesOnlyArray{i}(indexValuesOfNonZero{i}) = Zq{i}(indexValuesOfNonZero{i});
%             zzeroedWindowedWrinklesOnlyArray{i}(indexValuesOfNonZero{i}) = Zq{i}(indexValuesOfNonZero{i})-zMinForPlots;
        end
    end
end

for i = framesWithWrinkles
%     for x = floor(((MaximumBx{i}-MinimumBx{i})/tablingvaluex{i})/xWindowSize)*xWindowSize+1:numel(Zq{1}(1,:))
%         for y = floor(((MaximumBy{i}-MinimumBy{i})/tablingvaluey{i})/yWindowSize)*yWindowSize+1:numel(Zq{1}(:,1))
            logicalWindowedWrinklesOnlyArray{i}(floor(((MaximumBy{i}-MinimumBy{i})/tablingvaluey{i})/yWindowSize)*yWindowSize+1:numel(Zq{1}(:,1)),floor(((MaximumBx{i}-MinimumBx{i})/tablingvaluex{i})/xWindowSize)*xWindowSize+1:numel(Zq{1}(1,:))) = nansum(nansum(wrinklesOnlyArray{i}(floor(((MaximumBy{i}-MinimumBy{i})/tablingvaluey{i})/yWindowSize)*yWindowSize+1:numel(Zq{1}(:,1)),floor(((MaximumBx{i}-MinimumBx{i})/tablingvaluex{i})/xWindowSize)*xWindowSize+1:numel(Zq{1}(1,:)))));
            indexValuesOfNonZero{i} = find(logicalWindowedWrinklesOnlyArray{i});
            windowedWrinklesOnlyArray{i}(indexValuesOfNonZero{i}) = Zq{i}(indexValuesOfNonZero{i});
            zzeroedWindowedWrinklesOnlyArray{i}(indexValuesOfNonZero{i}) = Zq{i}(indexValuesOfNonZero{i})-zMinForPlots;
%         end
%     end
end

%%

% Plots the method for determining the maximum and minimum peak.
% Plots the original data set and a filtered version of that data set.
% Plots the detrended data set and a filtered version of the detrended data
% set.
% Plots the positive and negative peaks of the detrended filtered data set.
% Plots the maximum and corresponding two minimum peaks on the original
% data set.

 for i = 1:EndRunFiles-StartRunFiles+1;
 
     for j = ceil(((LastYVal{i}-FirstYVal{i})/2-MinimumBy{i})/tablingvaluey{i})%ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((floor((MaximumBy{i}-MinimumBy{i})/(tablingvaluey{i}))-1)/(3-1)):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})
         maximumWrinkleSurfacePlot{i} = figure('units','normalized','outerposition',[0 0 1 1]);
 %         maximumWrinkleSurfacePlot{i} = figure('Position',[0,0,1920,1080],'visible','off');
 
         h1 = subplot_tight(2,2,1,0.05);
         plot(ConstElementClippedXZAlongConstYDistance{i,j}(:,1),ConstElementClippedXZAlongConstYDistance{i,j}(:,3),'color',[135 206 250]/255);
         hold on
         plot(nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(:,1),nonDetrendedSmoothedConstElementClippedXZAlongConstYDistanceN{i,j}(:,2),'LineWidth',3);
         title({[names{i-1+StartRunFiles}],[(j-1)*tablingvaluey{i}+MinimumBy{i}]},'Interpreter','none');
         hold off
 
         h2 = subplot_tight(2,2,2,0.05);
         plot(DetrendedXZAlongConstYDistance{i,j}(:,1),DetrendedXZAlongConstYDistance{i,j}(:,3),'color',[135 206 250]/255);
         hold on
         plot(SmoothedConstElementClippedXZAlongConstYDistance{i,j}(:,1),SmoothedConstElementClippedXZAlongConstYDistance{i,j}(:,2),'LineWidth',3);
         plot(locs_MaxPeakXZConstY{i,j},pks_MaxPeakXZConstY{i,j},'wv','MarkerFaceColor','r');
         plot(locs_MinPeakXZConstY{i,j},-pks_MinPeakXZConstY{i,j},'wv','MarkerFaceColor','b');
     %     plot(sortedWrinkleArray{i,j}(end,5),sortedWrinkleArray{i,j}(end,7)-zMinForPlots,'wv','MarkerFaceColor','g');
         plot(sortedWrinkleArray{i,j}(end,5),SmoothedConstElementClippedXZAlongConstYDistance{i,j}(SmoothedConstElementClippedXZAlongConstYDistance{i,j}(:,1)==sortedWrinkleArray{i,j}(end,5),2),'wv','MarkerFaceColor',[0 0.75 0]);
         title([i;j]);
         hold off
 
         h3 = subplot_tight(2,2,3,0.05);
 
         plotxmin = -inf;
         plotxmax = inf;
         plotymin = -inf;
         plotymax = inf;
         plotzmin = 0;
         ceil(max(max(Zq{i}))-zMinForPlots)
         plotzmax = 2; % Set this equal to the print out value above
         plotzmax2 = 2;
         plotzmincolor = 0;
         plotzmaxcolor = 2;
         chosenj = j;
 
         % Plot flat and adjusted x based on initial position
         surf(Xg{i},Yg{i},Zq{i}-zMinForPlots)
         title(char(names(i-1+StartRunFiles)),'interpreter','none')
         xlabel('x [mm]')
         ylabel('y [mm]')
         zlabel('z [mm]')
         axis([plotxmin,plotxmax,plotymin,plotymax,plotzmin,plotzmax])
         caxis([plotzmincolor plotzmaxcolor])
         daspect([10 10 2])
         pbaspect([1 1 1])
         view([-5,30])
     %     view([0.5,90])
         c= colorbar;
         c.Label.String = 'z (mm)';
         shading interp
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%% Divide the lengths by the number of lines needed
         tempsy1{i} = size(Yg{i});
         tempsx1{i} = size(Xg{i});
         sxsize1{i} = floor(tempsx1{i}(1,2)/18); % 20 lines 
         sysize1{i} = floor(tempsy1{i}(1,1)/10); % 10 partitions
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         hold on
 
         % Make sure the dimensions are always correct based on surface plot data 
         % PLOTTING Lines in the X-Z plane
         for b1 = 1:sysize1{i}:tempsy1{i}(1,1)
             Y11s{i}(1:tempsy1{i}(1,2)) = Yg{i}(b1, 1);
             X11s{i}=Xg{i}(1,:);
             Z11s{i}=Zq{i}(b1,:)-zMinForPlots; 
             line(X11s{i},Y11s{i},Z11s{i},'color','black') ;
         end
         %     PLOTTING Lines in the Y-Z plane
         for c1 = 1:sxsize1{i}:tempsx1{i}(1,2)
             X21s{i}(1:tempsx1{i}(1,1))=Xg{i}(1, c1);
             Y21s{i}=Yg{i}(:, 1);
             Z21s{i}=Zq{i}(:, c1)-zMinForPlots;
             line(X21s{i},Y21s{i},Z21s{i},'color','black') ;
         end
         plot3(Xg{i}(chosenj,1:end),Yg{i}(chosenj,1:end),Zq{i}(chosenj,1:end)-zMinForPlots,'b','LineWidth',10);
         scatter3(wrinkleArrayCondensed{i}(:,1),wrinkleArrayCondensed{i}(:,2),wrinkleArrayCondensed{i}(:,3)-zMinForPlots,repmat(50,numel(wrinkleArrayCondensed{i}(:,1)),1),'filled',...
             'MarkerEdgeColor','r',...
             'MarkerFaceColor','r');
         scatter3(sortedWrinkleArray{i,j}(end,5),sortedWrinkleArray{i,j}(end,6),sortedWrinkleArray{i,j}(end,7)-zMinForPlots,repmat(150,numel(sortedWrinkleArray{i,j}(end,5)),1),'filled',...
             'MarkerEdgeColor',[0 0.75 0],...
             'MarkerFaceColor',[0 0.75 0]);
         hold off
 
         h4 = subplot_tight(2,2,4,0.05);
         plotxmin = -inf;
         plotxmax = inf;
         plotymin = -inf;
         plotymax = inf;
         plotzmin = 0;
         ceil(max(max(Zq{i}))-zMinForPlots)
         plotzmax = 2; % Set this equal to the print out value above
         plotzmax2 = 2;
         plotzmincolor = 0;
         plotzmaxcolor = 2;
         chosenj = j;
 
         % Plot flat and adjusted x based on initial position
         surf(Xg{i},Yg{i},Zq{i}-zMinForPlots)
         title(char(names(i-1+StartRunFiles)),'interpreter','none')
         xlabel('x [mm]')
         ylabel('y [mm]')
         zlabel('z [mm]')
         axis([plotxmin,plotxmax,plotymin,plotymax,plotzmin,plotzmax])
         caxis([plotzmincolor plotzmaxcolor])
         daspect([10 10 2])
         pbaspect([1 1 1])
         view([0.5,90])
     %     view([0.5,90])
         c= colorbar;
         c.Label.String = 'z (mm)';
         shading interp
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%% Divide the lengths by the number of lines needed
         tempsy1{i} = size(Yg{i});
         tempsx1{i} = size(Xg{i});
         sxsize1{i} = floor(tempsx1{i}(1,2)/18); % 20 lines 
         sysize1{i} = floor(tempsy1{i}(1,1)/10); % 10 partitions
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         hold on
 
         % Make sure the dimensions are always correct based on surface plot data 
         % PLOTTING Lines in the X-Z plane
         for b1 = 1:sysize1{i}:tempsy1{i}(1,1)
             Y11s{i}(1:tempsy1{i}(1,2)) = Yg{i}(b1, 1);
             X11s{i}=Xg{i}(1,:);
             Z11s{i}=Zq{i}(b1,:)-zMinForPlots; 
             line(X11s{i},Y11s{i},Z11s{i},'color','black') ;
         end
         %     PLOTTING Lines in the Y-Z plane
         for c1 = 1:sxsize1{i}:tempsx1{i}(1,2)
             X21s{i}(1:tempsx1{i}(1,1))=Xg{i}(1, c1);
             Y21s{i}=Yg{i}(:, 1);
             Z21s{i}=Zq{i}(:, c1)-zMinForPlots;
             line(X21s{i},Y21s{i},Z21s{i},'color','black') ;
         end
         surf(Xg{i},Yg{i},wrinklesOnlyArray{i}-zMinForPlots)
         plot3(Xg{i}(chosenj,1:end),Yg{i}(chosenj,1:end),Zq{i}(chosenj,1:end)-zMinForPlots,'b','LineWidth',5);
         scatter3(wrinkleArrayCondensed{i}(:,1),wrinkleArrayCondensed{i}(:,2),wrinkleArrayCondensed{i}(:,3)-zMinForPlots,repmat(50,numel(wrinkleArrayCondensed{i}(:,1)),1),'filled',...
             'MarkerEdgeColor','r',...
             'MarkerFaceColor','r');
         scatter3(sortedWrinkleArray{i,j}(end,5),sortedWrinkleArray{i,j}(end,6),sortedWrinkleArray{i,j}(end,7)-zMinForPlots,repmat(150,numel(sortedWrinkleArray{i,j}(end,5)),1),'filled',...
             'MarkerEdgeColor',[0 0.75 0],...
             'MarkerFaceColor',[0 0.75 0]);
         
         % Bounding box lines
         plot3(Xg{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Yg{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Zq{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end)-zMinForPlots,'b','LineWidth',5);
         plot3(Xg{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Yg{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Zq{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end)-zMinForPlots,'b','LineWidth',5);
         
         plot3(firstXArray{i}(:,1),firstXArray{i}(:,2),firstXArray{i}(:,3)-min(Zq{i}(:)),'b','LineWidth',5);
         plot3(lastXArray{i}(:,1),lastXArray{i}(:,2),lastXArray{i}(:,3)-min(Zq{i}(:)),'b','LineWidth',5);
 
         hold off
         
 %         saveas(maximumWrinkleSurfacePlot{i},[RootFileOutput '02_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_MatLab_Wrinkles_Method.png'],'png')
 
     end
 
 end

%%

% This plot shows an early version of the grid method

% figure;
%         plotxmin = -inf;
%         plotxmax = inf;
%         plotymin = -inf;
%         plotymax = inf;
%         plotzmin = 0;
%         ceil(max(max(Zq{i}))-zMinForPlots)
%         plotzmax = 2; % Set this equal to the print out value above
%         plotzmax2 = 2;
%         plotzmincolor = 0;
%         plotzmaxcolor = 2;
%         chosenj = j;
% 
%         % Plot flat and adjusted x based on initial position
%         sobject = surf(Xg{i},Yg{i},Zq{i}-zMinForPlots)
%         title(char(names(i-1+StartRunFiles)),'interpreter','none')
%         xlabel('x [mm]')
%         ylabel('y [mm]')
%         zlabel('z [mm]')
%         axis([plotxmin,plotxmax,plotymin,plotymax,plotzmin,plotzmax])
%         caxis([plotzmincolor plotzmaxcolor])
%         daspect([10 10 2])
%         pbaspect([1 1 1])
%         view([0.5,90])
%     %     view([0.5,90])
%         c= colorbar;
%         c.Label.String = 'z (mm)';
%         shading interp
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%% Divide the lengths by the number of lines needed
%         tempsy1{i} = size(Yg{i});
%         tempsx1{i} = size(Xg{i});
%         sxsize1{i} = floor(tempsx1{i}(1,2)/18); % 20 lines 
%         sysize1{i} = floor(tempsy1{i}(1,1)/10); % 10 partitions
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         hold on
% 
%         % Make sure the dimensions are always correct based on surface plot data 
%         % PLOTTING Lines in the X-Z plane
%         for b1 = 1:sysize1{i}:tempsy1{i}(1,1)
%             Y11s{i}(1:tempsy1{i}(1,2)) = Yg{i}(b1, 1);
%             X11s{i}=Xg{i}(1,:);
%             Z11s{i}=Zq{i}(b1,:)-zMinForPlots; 
%             line(X11s{i},Y11s{i},Z11s{i},'color','black') ;
%         end
%         %     PLOTTING Lines in the Y-Z plane
%         for c1 = 1:sxsize1{i}:tempsx1{i}(1,2)
%             X21s{i}(1:tempsx1{i}(1,1))=Xg{i}(1, c1);
%             Y21s{i}=Yg{i}(:, 1);
%             Z21s{i}=Zq{i}(:, c1)-zMinForPlots;
%             line(X21s{i},Y21s{i},Z21s{i},'color','black') ;
%         end
%         surf(Xg{i},Yg{i},wrinklesOnlyArray{i}-zMinForPlots)
%         surf(Xg{i},Yg{i},zzeroedWindowedWrinklesOnlyArray{i})
% %         plot3(Xg{i}(chosenj,1:end),Yg{i}(chosenj,1:end),Zq{i}(chosenj,1:end)-zMinForPlots,'b','LineWidth',10);
% %         scatter3(wrinkleArrayCondensed{i}(:,1),wrinkleArrayCondensed{i}(:,2),wrinkleArrayCondensed{i}(:,3)-zMinForPlots,repmat(100,numel(wrinkleArrayCondensed{i}(:,1)),1),'filled',...
% %             'MarkerEdgeColor','r',...
% %             'MarkerFaceColor','r');
% %         scatter3(sortedWrinkleArray{i,j}(end,5),sortedWrinkleArray{i,j}(end,6),sortedWrinkleArray{i,j}(end,7)-zMinForPlots,repmat(300,numel(sortedWrinkleArray{i,j}(end,5)),1),'filled',...
% %             'MarkerEdgeColor',[0 0.75 0],...
% %             'MarkerFaceColor',[0 0.75 0]);
%         hold off

%%

clearvars centroidValuesMirroredExpanded centroidValuesMirrored centroidPixelValues centroidValues

% http://blogs.mathworks.com/steve/2009/02/27/using-ismember-with-the-output-of-regionprops/

% centroidValuesMirroredExpanded takes the logical array with the wrinkle
% positions and applies the region props function which returns a 'struct'
% with the centroid information.

% centroidValuesMirrored takes all the centroid values and puts them into a
% single array.

% centroidPixelValues this is a hold over, I used to think that the values
% were mirrored about the y axis, but it doesn't look like they are so just
% ignore this array.

% centroidValues converts the pixel values returned from the computation
% and converts them to x and y coordinates based on the tabling values.

for i = framesWithWrinkles
    centroidValuesMirroredExpanded{i} = regionprops(logicalWindowedWrinklesOnlyArray{i},'centroid');
    centroidValuesMirrored{i} = cat(1, centroidValuesMirroredExpanded{i}.Centroid);
%     centroidPixelValues{i} = [centroidValuesMirrored{i}(:,1) numel(logicalWindowedWrinklesOnlyArray{i}(:,1))-centroidValuesMirrored{i}(:,2)];
    centroidPixelValues{i} = [centroidValuesMirrored{i}(:,1) centroidValuesMirrored{i}(:,2)];
    centroidValues{i} = [tablingvaluex{i}*centroidPixelValues{i}(:,1) tablingvaluey{i}*centroidPixelValues{i}(:,2)];
end

%%

% wrinkleLabels is the identifier for each of the zones.  The length of
% which is the number of wrinkles.  This is an array with each point being
% equal to the closest wrinkle denoted by a different value.  A map of 0's
% 1's for points closest to 1, 2's for points closest to wrinkle 2, etc.

for i = framesWithWrinkles
    wrinkleLabels{i} = bwlabel(logicalWindowedWrinklesOnlyArray{i});
end

%%

% Takes the individual wrinkle label then finds the match between that
% label and the map generated above, this is the x&y position of the
% containted wrinkle.  This will be continued below.

for i = framesWithWrinkles
    for k = 1:length(centroidValues{i})
        wrinkleZoneX{i,k} = Xg{i}(find(ismember(wrinkleLabels{i}, k)));
        wrinkleZoneY{i,k} = Yg{i}(find(ismember(wrinkleLabels{i}, k)));
    end
end

%%

% New method.  Old method used the max z value in the wrinkle zone.  This
% would give an error because the max z value did not necessarily
% correlate to a maximum wrinkle peak.  IE. areas with higher curvature
% would have a high z-value but the minimum peaks would also be high
% resulting in a small wrinkle relative to the substrate.

% New method finds all of the maximum peaks (calculated above, hence larger
% than the prominence value) which are contained in the wrinkleLabels
% region, k.  Take the maximum of those wrinkle heights, (along with the
% y-position) to determine the wrinkleLabels{k} maximum wrinkle height and
% width.

clear containedWrinklesArray

% Needed an expansion factor because some of the maximum wrinkle heights in
% x were 1 discrete step larger than the window resulting in no wrinkles
% being assigned to a wrinkleLabel with small boundaries.

expansionsFactorForContainedWrinkles = tablingvaluex{1};

% The logical applied here just takes the x&y bounds and finds the wrinkle
% from the condensed set which satisfy the x&y constraints.  IE the peaks
% inside that boundary.

for i = framesWithWrinkles
        for k = 1:length(centroidValues{i}(:,1))
                containedWrinklesArray{i,k} = wrinkleArrayCondensed{i}(wrinkleArrayCondensed{i}(:,1)>=min(wrinkleZoneX{i,k})-expansionsFactorForContainedWrinkles & wrinkleArrayCondensed{i}(:,1)<=max(wrinkleZoneX{i,k})+expansionsFactorForContainedWrinkles & wrinkleArrayCondensed{i}(:,2)>=min(wrinkleZoneY{i,k})-expansionsFactorForContainedWrinkles & wrinkleArrayCondensed{i}(:,2)<=max(wrinkleZoneY{i,k})+expansionsFactorForContainedWrinkles,:);
        end
end

%%

% Now takes the maximum wrinkle height from the containedWrinklesArray for
% each frame and wrinkleLabel k and finds the corresponding x,y,z value.

clear maxWrinkleHeightVal maxWrinkleHeightInd maxWrinkleXVal maxWrinkleYVal maxWrinkleZVal wrinkleMetricsUnbounded

for i = framesWithWrinkles
    wrinkleMetricsUnbounded{i} = [];
    for k = 1:length(centroidValues{i}(:,1))
        [maxWrinkleHeightVal{i,k} maxWrinkleHeightInd{i,k}] = max(containedWrinklesArray{i,k}(:,4));
        maxWrinkleXVal{i,k} = containedWrinklesArray{i,k}(maxWrinkleHeightInd{i,k},1);
        maxWrinkleYVal{i,k} = containedWrinklesArray{i,k}(maxWrinkleHeightInd{i,k},2);
        maxWrinkleZVal{i,k} = containedWrinklesArray{i,k}(maxWrinkleHeightInd{i,k},3);
        wrinkleMetricsUnbounded{i} = [wrinkleMetricsUnbounded{i}; k containedWrinklesArray{i,k}(maxWrinkleHeightInd{i,k},4:5) containedWrinklesArray{i,k}(maxWrinkleHeightInd{i,k},8) containedWrinklesArray{i,k}(maxWrinkleHeightInd{i,k},6:7) containedWrinklesArray{i,k}(maxWrinkleHeightInd{i,k},1:3)];
    end
end

%%

clear xMinBoundary xMaxBoundary yMinBoundary yMaxBoundary wrinkleMetrics

% Work around to remove some of the wrinkles at the boundaries causing
% issues because of the insufficient length to calculate a discontinuity
% between the polycarbonate and the prepreg.

% Set these as small as possible without resulting in an error to include
% as many wrinkles as the maxProminence allows.

xMinBoundary = -inf;
xMaxBoundary = inf;
yMinBoundary = -inf;
yMaxBoundary = inf;

for i = framesWithWrinkles
    wrinkleMetrics{i} = wrinkleMetricsUnbounded{i}(wrinkleMetricsUnbounded{i}(:,6) >= xMinBoundary & wrinkleMetricsUnbounded{i}(:,6) <= xMaxBoundary & wrinkleMetricsUnbounded{i}(:,7) >= yMinBoundary & wrinkleMetricsUnbounded{i}(:,7) <= yMaxBoundary,:);
end

%%

clear colorTableSorted sortingColorValue sortingColorIndex

for i = 1:EndRunFiles-StartRunFiles+1;

colorTableSorted{i} = autumn(1000);%length(wrinkleMetrics{i}));

    for k = 1:length(wrinkleMetrics{i}(:,1));
        [cColour{i,k} cIndex{i,k}] = min(abs((1-colorTableSorted{i}(:,2))-wrinkleMetrics{i}(k,2)/0.55));
        wColour{i,k} = colorTableSorted{i}(cIndex{i,k},:);
    end

% [cColour cIndex] = min(abs(colorTableSorted{i}(:,2)-wrinkleMetrics{1}(1,2)/0.4))
% [sortingColorValue{i},sortingColorIndex{i}] = sort(wrinkleMetrics{i}(:,2),'descend');
% colorTable{i} = colorTableSorted{i};%(sortingColorIndex{i},:);

    meanFirstX{i} = mean([FirstX{i,ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})}]);
    meanLastX{i} = mean([LastX{i,ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}):floor((LastYVal{i}-MinimumBy{i})/tablingvaluey{i})}]);

    zShiftPlot = 0.3;
    polycarbonateZAverage{i} = mean(Zq{i}(((Yg{i}>LastYVal{i}) | (Yg{i}<FirstYVal{i}) | (Xg{i}<meanFirstX{i}) | (Xg{i}>meanLastX{i}))))-zShiftPlot;

    plotzmin = 0;
    ceil(max(max(Zq{i}))-min(min(Zq{i})))
    plotzmax = 1; % Set this equal to the print out value above
    plotzmax2 = 1;

    plotzmincolor = 0;
    plotzmaxcolor = 1;
    
%     wrinkleLocationIdentificationFigure{i} = figure('units','normalized','outerposition',[0 0 1 1]);
    wrinkleLocationIdentificationFigure{i} = figure('Position',[0,0,1920,1080],'visible','off');
    %     h1 = subplot(2,3,[1,2,4,5]);
        surf(Xg{i},Yg{i},Zq{i}-polycarbonateZAverage{i})
    %     title(char(names(j-1+StartRunFiles)),'interpreter','none')
        xlabel('x [mm]')
        ylabel('y [mm]')
        ylabh = get(gca,'ylabel');
        set(ylabh,'Units','normalized');
        set(ylabh,'position',get(ylabh,'position') - [0.35 -0.3 0]);
        xlabh = get(gca,'xlabel');
        set(xlabh,'Units','normalized');
        set(xlabh,'position',get(xlabh,'position') - [0.1 0.1 0]);
    %     zlabel('z [mm]')
    %     axis([plotxmin,plotxmax,plotymin,plotymax,plotzmin,plotzmax])
%         colormap(flipud(autumn));
        axis([0,160,0,80,plotzmin,plotzmax])
%         axis([0,160,0,80,0,0.6])
        xticks(0:20:160)
        yticks(0:20:80)
        set(gca,'ztick',[])
    %     xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi])
        caxis([plotzmincolor plotzmaxcolor])
%         caxis([0 0.6])
        daspect([10 10 2])
        pbaspect([1 1 1])
        view([-10,30])
%         view([0.5,90])
        c= colorbar;
%         c= colorbar('Ticks',[0,0.3,0.6]);
        c.Label.String = 'z [mm]';
        set(c,'position',[.9 .45 .0175 .3])
        shading interp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Divide the lengths by the number of lines needed
        tempsy1{i} = size(Yg{i});
        tempsx1{i} = size(Xg{i});
        sxsize1{i} = round(160/157*3141/16);%floor(tempsx1{j}(1,2)/18); % 20 lines 
        sysize1{i} = round(80/80*1601/8);%floor(tempsy1{j}(1,1)/10); % 10 partitions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on

        % Make sure the dimensions are always correct based on surface plot data 
        % PLOTTING Lines in the X-Z plane
        for b1 = 1:sysize1{i}:tempsy1{i}(1,1)
            Y11s{i}(1:tempsy1{i}(1,2)) = Yg{i}(b1, 1);
            X11s{i}=Xg{i}(1,:);
            Z11s{i}=Zq{i}(b1,:)-polycarbonateZAverage{i}; 
            line(X11s{i},Y11s{i},Z11s{i},'color','black') ;
        end
        %     PLOTTING Lines in the Y-Z plane
        for c1 = 1:sxsize1{i}:tempsx1{i}(1,2)
            X21s{i}(1:tempsx1{i}(1,1))=Xg{i}(1, c1);
            Y21s{i}=Yg{i}(:, 1);
            Z21s{i}=Zq{i}(:, c1)-polycarbonateZAverage{i};
            line(X21s{i},Y21s{i},Z21s{i},'color','black') ;
        end
%         surf(Xg{i},Yg{i},wrinklesOnlyArray{i}-polycarbonateZAverage{i})
%         surf(Xg{i},Yg{i},zzeroedWindowedWrinklesOnlyArray{i}-zShiftPlot)
        
        for k = 1:length(wrinkleMetrics{i}(:,1));
        scatter3(Xg{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1)))),Yg{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1)))),Zq{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1))))-polycarbonateZAverage{i},5,'filled',...
            'MarkerEdgeColor',wColour{i,k},...
            'MarkerFaceColor',wColour{i,k});
        end
        
%         for k = 1:length(wrinkleMetrics{i});
%         scatter3(Xg{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1)))),Yg{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1)))),Zq{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1))))-min(min(Zq{i})),5,'filled',...
%             'MarkerEdgeColor',colorTable{i}(ismember(wrinkleMetrics{i}(sortingColorIndex{i},1),wrinkleMetrics{i}(k,1)),:),...
%             'MarkerFaceColor',colorTable{i}(ismember(wrinkleMetrics{i}(sortingColorIndex{i},1),wrinkleMetrics{i}(k,1)),:));
%         end
        
        % Wrinkle x-y gridlines
        greenLineFactor = 3;
        for m = 1:length(wrinkleMetrics{i}(:,1));
            plot3(Xg{i}(floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}-greenLineFactor)/tablingvaluey{i}):floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}+greenLineFactor)/tablingvaluey{i}),ceil((wrinkleMetrics{i}(m,7)-MinimumBy{i})/tablingvaluey{i})),Yg{i}(floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}-greenLineFactor)/tablingvaluey{i}):floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}+greenLineFactor)/tablingvaluey{i}),ceil((wrinkleMetrics{i}(m,7)-MinimumBy{i})/tablingvaluey{i})),Zq{i}(floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}-greenLineFactor)/tablingvaluey{i}):floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}+greenLineFactor)/tablingvaluey{i}),ceil((wrinkleMetrics{i}(m,7)-MinimumBy{i})/tablingvaluey{i}))-polycarbonateZAverage{i},'g','LineWidth',1);
            plot3(Xg{i}(ceil((wrinkleMetrics{i}(m,8)-MinimumBy{i})/tablingvaluex{i}),floor((wrinkleMetrics{i}(m,7)+MinimumBx{i}-greenLineFactor)/tablingvaluex{i}):ceil((wrinkleMetrics{i}(m,7)+MinimumBx{i}+greenLineFactor)/tablingvaluex{i})),Yg{i}(ceil((wrinkleMetrics{i}(m,8)-MinimumBy{i})/tablingvaluex{i}),floor((wrinkleMetrics{i}(m,7)+MinimumBx{i}-greenLineFactor)/tablingvaluex{i}):ceil((wrinkleMetrics{i}(m,7)+MinimumBx{i}+greenLineFactor)/tablingvaluex{i})),Zq{i}(ceil((wrinkleMetrics{i}(m,8)-MinimumBy{i})/tablingvaluex{i}),floor((wrinkleMetrics{i}(m,7)+MinimumBx{i}-greenLineFactor)/tablingvaluex{i}):ceil((wrinkleMetrics{i}(m,7)+MinimumBx{i}+greenLineFactor)/tablingvaluex{i}))-polycarbonateZAverage{i},'g','LineWidth',1);
        end

%         % Wrinkle text identifier
%         for m = 1:length(wrinkleMetrics{i});
% %             t{m} = text(wrinkleMetrics{i}(m,7)-3,wrinkleMetrics{i}(m,8)-1,wrinkleMetrics{i}(m,9)-min(Zq{i}(:))+0.2,num2str(wrinkleMetrics{i}(m,1)));
%             t{m} = text(wrinkleMetrics{i}(m,7)-1,wrinkleMetrics{i}(m,8)-1,wrinkleMetrics{i}(m,9)-min(Zq{i}(:))+0.2,num2str(wrinkleMetrics{i}(m,1)));
%             t{m}.FontSize = 7;
%             t{m}.Margin = 1;
%             t{m}.BackgroundColor = 'w';
%             t{m}.EdgeColor = 'k';
%             t{m}.HorizontalAlignment = 'center';
%             t{m}.VerticalAlignment = 'middle';
%         end

%         % Bounding box lines
%         plot3(Xg{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Yg{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Zq{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end)-min(min(Zq{i})),'b','LineWidth',5);
%         plot3(Xg{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Yg{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Zq{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end)-min(min(Zq{i})),'b','LineWidth',5);
%         
%         plot3(firstXArray{i}(:,1),firstXArray{i}(:,2),firstXArray{i}(:,3)-min(Zq{i}(:)),'b','LineWidth',5);
%         plot3(lastXArray{i}(:,1),lastXArray{i}(:,2),lastXArray{i}(:,3)-min(Zq{i}(:)),'b','LineWidth',5);

        hold off
        set(gca, 'FontSize', 14)
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 22 8]); %x_width=30cm y_width=25cm
%         set(gcf, 'PaperPosition', [0 0 44 16]); %x_width=30cm y_width=25cm
%         set(gcf, 'OuterPosition',[0 0 0.1 0.1]);
        
%         saveas(wrinkleLocationIdentificationFigure{i},[RootFileOutput '05_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_Wrinkle_Identifier.png'],'png')
%         print(wrinkleLocationIdentificationFigure{i},[RootFileOutput '05_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_Wrinkle_Identifier.png'],'-dpng','-r300')
%         print(wrinkleLocationIdentificationFigure{i},[RootFileOutput '06_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_Wrinkle_Identifier_Full.png'],'-dpng','-r300')
        print(wrinkleLocationIdentificationFigure{i},[RootFileOutput '07_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_Wrinkle_Full.png'],'-dpng','-r300')

end

%%

clear colorTableSorted sortingColorValue sortingColorIndex

for i = 1%framesWithWrinkles;

colorTableSorted{i} = autumn(1000);%length(wrinkleMetrics{i}));

    for k = 1:length(wrinkleMetrics{i}(:,1));
        [cColour{i,k} cIndex{i,k}] = min(abs((1-colorTableSorted{i}(:,2))-wrinkleMetrics{i}(k,2)/0.55));
        wColour{i,k} = colorTableSorted{i}(cIndex{i,k},:);
    end
    
% [cColour cIndex] = min(abs(colorTableSorted{i}(:,2)-wrinkleMetrics{1}(1,2)/0.4))
% [sortingColorValue{i},sortingColorIndex{i}] = sort(wrinkleMetrics{i}(:,2),'descend');
% colorTable{i} = colorTableSorted{i};%(sortingColorIndex{i},:);

    plotzmin = 0;
    ceil(max(max(Zq{i}))-min(min(Zq{i})))
    plotzmax = 1; % Set this equal to the print out value above
    plotzmax2 = 1;

    plotzmincolor = 0;
    plotzmaxcolor = 1;
    
%     wrinkleLocationIdentificationFigure{i} = figure('units','normalized','outerposition',[0 0 1 1]);
    wrinkleLocationIdentificationFigure{i} = figure('Position',[0,0,1920,1080],'visible','off');
    %     h1 = subplot(2,3,[1,2,4,5]);
        surf(Xg{i},Yg{i},Zq{i}-polycarbonateZAverage{i})
    %     title(char(names(j-1+StartRunFiles)),'interpreter','none')
        xlabel('x [mm]')
        ylabel('y [mm]')
%         ylabh = get(gca,'ylabel');
%         set(ylabh,'Units','normalized');
%         set(ylabh,'position',get(ylabh,'position') - [0.35 -0.3 0]);
%         xlabh = get(gca,'xlabel');
%         set(xlabh,'Units','normalized');
%         set(xlabh,'position',get(xlabh,'position') - [0.1 0.1 0]);
    %     zlabel('z [mm]')
    %     axis([plotxmin,plotxmax,plotymin,plotymax,plotzmin,plotzmax])
%         colormap(flipud(autumn));
        axis([0,160,0,80,plotzmin,plotzmax])
%         axis([0,160,0,80,0,0.6])
        xticks(0:20:160)
        yticks(0:20:80)
        set(gca,'ztick',[])
    %     xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi])
        caxis([plotzmincolor plotzmaxcolor])
%         caxis([0 0.6])
        daspect([10 10 2])
        pbaspect([1 1 1])
%         view([-10,30])
        view([0.5,90])
        c= colorbar;
%         c= colorbar('Ticks',[0,0.3,0.6]);
        c.Label.String = 'z [mm]';
        set(c,'position',[.9 .45 .0175 .3])
        shading interp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Divide the lengths by the number of lines needed
        tempsy1{i} = size(Yg{i});
        tempsx1{i} = size(Xg{i});
        sxsize1{i} = round(160/157*3141/16);%floor(tempsx1{j}(1,2)/18); % 20 lines 
        sysize1{i} = round(80/80*1601/8);%floor(tempsy1{j}(1,1)/10); % 10 partitions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on

        % Make sure the dimensions are always correct based on surface plot data 
        % PLOTTING Lines in the X-Z plane
        for b1 = 1:sysize1{i}:tempsy1{i}(1,1)
            Y11s{i}(1:tempsy1{i}(1,2)) = Yg{i}(b1, 1);
            X11s{i}=Xg{i}(1,:);
            Z11s{i}=Zq{i}(b1,:)-polycarbonateZAverage{i}; 
            line(X11s{i},Y11s{i},Z11s{i},'color','black') ;
        end
        %     PLOTTING Lines in the Y-Z plane
        for c1 = 1:sxsize1{i}:tempsx1{i}(1,2)
            X21s{i}(1:tempsx1{i}(1,1))=Xg{i}(1, c1);
            Y21s{i}=Yg{i}(:, 1);
            Z21s{i}=Zq{i}(:, c1)-polycarbonateZAverage{i};
            line(X21s{i},Y21s{i},Z21s{i},'color','black') ;
        end
%         surf(Xg{i},Yg{i},wrinklesOnlyArray{i}-min(min(Zq{i})))
%         surf(Xg{i},Yg{i},zzeroedWindowedWrinklesOnlyArray{i})
        
        for k = 1:length(wrinkleMetrics{i}(:,1));
        scatter3(Xg{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1)))),Yg{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1)))),Zq{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1))))-polycarbonateZAverage{i},5,'filled',...
            'MarkerEdgeColor',wColour{i,k},...
            'MarkerFaceColor',wColour{i,k});
        end
        
%         for k = 1:length(wrinkleMetrics{i});
%         scatter3(Xg{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1)))),Yg{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1)))),Zq{i}(find(ismember(wrinkleLabels{i}, wrinkleMetrics{i}(k,1))))-min(min(Zq{i})),5,'filled',...
%             'MarkerEdgeColor',colorTable{i}(ismember(wrinkleMetrics{i}(sortingColorIndex{i},1),wrinkleMetrics{i}(k,1)),:),...
%             'MarkerFaceColor',colorTable{i}(ismember(wrinkleMetrics{i}(sortingColorIndex{i},1),wrinkleMetrics{i}(k,1)),:));
%         end
        
        % Wrinkle x-y gridlines
        greenLineFactor = 3;
        for m = 1:length(wrinkleMetrics{i}(:,1));
            plot3(Xg{i}(floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}-greenLineFactor)/tablingvaluey{i}):floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}+greenLineFactor)/tablingvaluey{i}),ceil((wrinkleMetrics{i}(m,7)-MinimumBy{i})/tablingvaluey{i})),Yg{i}(floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}-greenLineFactor)/tablingvaluey{i}):floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}+greenLineFactor)/tablingvaluey{i}),ceil((wrinkleMetrics{i}(m,7)-MinimumBy{i})/tablingvaluey{i})),Zq{i}(floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}-greenLineFactor)/tablingvaluey{i}):floor((wrinkleMetrics{i}(m,8)+MinimumBy{i}+greenLineFactor)/tablingvaluey{i}),ceil((wrinkleMetrics{i}(m,7)-MinimumBy{i})/tablingvaluey{i}))-polycarbonateZAverage{i},'g','LineWidth',1);
            plot3(Xg{i}(ceil((wrinkleMetrics{i}(m,8)-MinimumBy{i})/tablingvaluex{i}),floor((wrinkleMetrics{i}(m,7)+MinimumBx{i}-greenLineFactor)/tablingvaluex{i}):ceil((wrinkleMetrics{i}(m,7)+MinimumBx{i}+greenLineFactor)/tablingvaluex{i})),Yg{i}(ceil((wrinkleMetrics{i}(m,8)-MinimumBy{i})/tablingvaluex{i}),floor((wrinkleMetrics{i}(m,7)+MinimumBx{i}-greenLineFactor)/tablingvaluex{i}):ceil((wrinkleMetrics{i}(m,7)+MinimumBx{i}+greenLineFactor)/tablingvaluex{i})),Zq{i}(ceil((wrinkleMetrics{i}(m,8)-MinimumBy{i})/tablingvaluex{i}),floor((wrinkleMetrics{i}(m,7)+MinimumBx{i}-greenLineFactor)/tablingvaluex{i}):ceil((wrinkleMetrics{i}(m,7)+MinimumBx{i}+greenLineFactor)/tablingvaluex{i}))-polycarbonateZAverage{i},'g','LineWidth',1);
        end

        % Wrinkle text identifier
        for m = 1:length(wrinkleMetrics{i}(:,1));
%             t{m} = text(wrinkleMetrics{i}(m,7)-3,wrinkleMetrics{i}(m,8)-1,wrinkleMetrics{i}(m,9)-min(Zq{i}(:))+0.2,num2str(wrinkleMetrics{i}(m,1)));
            t{m} = text(wrinkleMetrics{i}(m,7)-3,wrinkleMetrics{i}(m,8)-3,wrinkleMetrics{i}(m,9)-min(Zq{i}(:))+0.2,num2str(wrinkleMetrics{i}(m,1)));
            t{m}.FontSize = 10;
            t{m}.Margin = 1;
            t{m}.BackgroundColor = 'w';
            t{m}.EdgeColor = 'k';
            t{m}.HorizontalAlignment = 'center';
            t{m}.VerticalAlignment = 'middle';
        end

%         % Bounding box lines
%         plot3(Xg{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Yg{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Zq{i}(ceil((FirstYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end)-min(min(Zq{i})),'b','LineWidth',5);
%         plot3(Xg{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Yg{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end),Zq{i}(ceil((LastYVal{i}-MinimumBy{i})/tablingvaluey{i}),1:end)-min(min(Zq{i})),'b','LineWidth',5);
%         
%         plot3(firstXArray{i}(:,1),firstXArray{i}(:,2),firstXArray{i}(:,3)-min(Zq{i}(:)),'b','LineWidth',5);
%         plot3(lastXArray{i}(:,1),lastXArray{i}(:,2),lastXArray{i}(:,3)-min(Zq{i}(:)),'b','LineWidth',5);

        hold off
        set(gca, 'FontSize', 14)
        set(gcf, 'PaperUnits', 'centimeters');
%         set(gcf, 'PaperPosition', [0 0 22 8]); %x_width=30cm y_width=25cm
        set(gcf, 'PaperPosition', [0 0 32 14]); %x_width=30cm y_width=25cm
%         set(gcf, 'OuterPosition',[0 0 0.1 0.1]);
        
%         saveas(wrinkleLocationIdentificationFigure{i},[RootFileOutput '05_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_Wrinkle_Identifier.png'],'png')
%         print(wrinkleLocationIdentificationFigure{i},[RootFileOutput '05_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_Wrinkle_Identifier.png'],'-dpng','-r300')
        print(wrinkleLocationIdentificationFigure{i},[RootFileOutput '06_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_Wrinkle_Identifier_Full.png'],'-dpng','-r300')
%         print(wrinkleLocationIdentificationFigure{i},[RootFileOutput '07_' strrep(char(names(i-1+StartRunFiles)),'.xyz','') '_Wrinkle_Full.png'],'-dpng','-r300')

end

%%

% Output the values to Excel files as required.

% for i = framesWithWrinkles%1:EndRunFiles-StartRunFiles+1
%     outputWrinkleAtFrameSheet{i} = cell2table(num2cell([wrinkleMetrics{i}(:,1) repmat(i+StartRunFiles-1,length(wrinkleMetrics{i}(:,1)),1) wrinkleMetrics{i}(:,2:end) repmat(meanCurvature{i},length(wrinkleMetrics{i}(:,1)),1) repmat(stdCurvature{i},length(wrinkleMetrics{i}(:,1)),1)]),...
%         'VariableNames',{'Wrinkle_ID','Frame','Height','Width','Path_length','x_pos_left', 'x_pos_right', 'x_pos_center', 'y_pos_center', 'z_pos_center', 'Mean_curvature', 'STD_curvature'});
%     writetable(outputWrinkleAtFrameSheet{i},[RootFileOutput '07_' ExperimentDateCode '_' ExperimentNumber '_MatLab_Wrinkle_Metrics_' names{StartRunFiles}(end-10:end-9) '_' names{EndRunFiles}(end-10:end-9) '.xlsx'],'Sheet',i,'Range','A1')
% end

%%

toc
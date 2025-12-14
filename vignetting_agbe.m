clear; close all; clc;


ieInit;   


%imgFile = '../images/image2.png';  
imageSize  = 512;   % total image size (pixels)
checkerSize = 32;   % size of each square (pixels)

[X, Y] = meshgrid(1:imageSize, 1:imageSize);
checker = mod(floor(X/checkerSize) + floor(Y/checkerSize), 2);

originalImage = checker * 100 + 100;   

originalImage = mat2gray(originalImage); 

rgbIn = repmat(originalImage, [1 1 3]);  
[H, W, ~] = size(rgbIn);
side = min(H, W);             

rowStart = floor((H - side)/2) + 1;
colStart = floor((W - side)/2) + 1;
rowEnd   = rowStart + side - 1;
colEnd   = colStart + side - 1;

rgbIn = rgbIn(rowStart:rowEnd, colStart:colEnd, :);   


maxDim = 512;
if side > maxDim
    scale = maxDim / side;
    rgbIn = imresize(rgbIn, scale);
end

scene = sceneFromFile(rgbIn, 'rgb', 100);  
scene = sceneSet(scene, 'hfov', 70);

oi = oiCreate('diffraction limited');
oi = oiSet(oi, 'optics fnumber', 2.8);
oi = oiCompute(oi, scene);

oiOriginal = oi;
oiOriginal = oiSet(oiOriginal, 'name', 'OI - No Vignetting');


photons = oiGet(oiOriginal, 'photons');
[rows, cols, nWave] = size(photons);

spatialSupport = oiGet(oiOriginal, 'spatial support', 'm');

if iscell(spatialSupport)
    xSupport = spatialSupport{1};
    ySupport = spatialSupport{2};
elseif isnumeric(spatialSupport)
    if ndims(spatialSupport) == 3 && size(spatialSupport,3) == 2
        xSupport = spatialSupport(:,:,1);
        ySupport = spatialSupport(:,:,2);
    else
        error('Unexpected numeric format for spatial support from oiGet.');
    end
elseif isstruct(spatialSupport)
    if isfield(spatialSupport, 'x') && isfield(spatialSupport, 'y')
        xSupport = spatialSupport.x;
        ySupport = spatialSupport.y;
    elseif isfield(spatialSupport, 'xsupport') && isfield(spatialSupport, 'ysupport')
        xSupport = spatialSupport.xsupport;
        ySupport = spatialSupport.ysupport;
    else
        error('Unknown struct format for spatial support from oiGet.');
    end
else
    error('Unsupported type returned by oiGet(''spatial support'').');
end

if ~isequal(size(xSupport,1), rows) || ~isequal(size(xSupport,2), cols)
    error('Size mismatch: spatial support and photons do not have the same spatial dimensions.');
end

r = sqrt(xSupport.^2 + ySupport.^2);

hfov_deg = sceneGet(scene, 'hfov');   
hfov_rad = deg2rad(hfov_deg);

centerRow = round(rows/2);
r_edge = max(abs(xSupport(centerRow, :)));  

f_est = r_edge / tan(hfov_rad/2);

theta = atan(r ./ f_est);

vignetteStrength = 0.75;
thetaEff = vignetteStrength * theta;

vMap2D = cos(thetaEff).^4;
vMap2D(vMap2D < 0) = 0;   

vMap3D = repmat(vMap2D, [1 1 nWave]);

photonsVig   = photons .* vMap3D;
oiVignetted  = oiSet(oiOriginal, 'photons', photonsVig);
oiVignetted  = oiSet(oiVignetted, 'name', 'OI - With Angular cos^4 Vignetting');

fov = sceneGet(scene, 'hfov');

sensor = sensorCreate('bayer (rggb)');
sensor = sensorSetSizeToFOV(sensor, fov, oiOriginal);

sensorOrig = sensorCompute(sensor, oiOriginal);
ipOrig     = ipCreate();
ipOrig     = ipCompute(ipOrig, sensorOrig);

sensorVig = sensorCompute(sensor, oiVignetted);
ipVig     = ipCreate();
ipVig     = ipCompute(ipVig, sensorVig);

rgbOrig = ipGet(ipOrig, 'srgb');   % original, no vignetting
rgbVig  = ipGet(ipVig,  'srgb');   % vignetted

rgbOrig = double(rgbOrig);
rgbVig  = double(rgbVig);

origGray = 0.2989 * rgbOrig(:,:,1) + 0.5870 * rgbOrig(:,:,2) + 0.1140 * rgbOrig(:,:,3);
vigGray  = 0.2989 * rgbVig(:,:,1)  + 0.5870 * rgbVig(:,:,2)  + 0.1140 * rgbVig(:,:,3);

imgNorm = vigGray / max(vigGray(:));   % normalize to [0,1]
[rowsIP, colsIP] = size(imgNorm);
vignetteMapEst = ones(rowsIP, colsIP); % initial flat estimate

nIter = 50;

for iter = 1:nIter
    tempDevig = imgNorm ./ (vignetteMapEst + 1e-6);
    tempDevig = min(max(tempDevig, 0), 1);
    
    gaussianSigma = max(rowsIP, colsIP) / 6;
    H = fspecial('gaussian', [51 51], gaussianSigma);
    localIntensity = imfilter(tempDevig, H, 'replicate');
    localIntensity(localIntensity < 0.01) = 0.01;
    
    vignetteMapEst = localIntensity / max(localIntensity(:));
end

vignetteMapEst(vignetteMapEst < 0.01) = 0.01;
devigNorm = imgNorm ./ vignetteMapEst;
devigNorm = min(max(devigNorm, 0), 1);

devigGray = devigNorm * max(vigGray(:));


origNorm = origGray / max(origGray(:));
vigNorm  = vigGray  / max(vigGray(:));
devigNormForMetrics = devigGray / max(devigGray(:));

mse     = mean((origNorm(:) - devigNormForMetrics(:)).^2);
psnrVal = 10 * log10(1^2 / mse);   

residual = devigNormForMetrics - origNorm;

hLap      = fspecial('laplacian', 0.2);
lapOrig   = imfilter(origNorm, hLap, 'replicate');
lapVig    = imfilter(vigNorm,  hLap, 'replicate');
lapDevig  = imfilter(devigNormForMetrics, hLap, 'replicate');

sharpOrig  = var(lapOrig(:));
sharpVig   = var(lapVig(:));
sharpDevig = var(lapDevig(:));

fprintf('=== De-Vignetting Metrics (AGBE, ISETCam, Angular cos^4) ===\n');
fprintf('PSNR (De-vignetted vs Original) : %.2f dB\n', psnrVal);
fprintf('Sharpness (variance of Laplacian):\n');
fprintf('  Original (no vignetting)  : %.4g\n', sharpOrig);
fprintf('  Vignetted                 : %.4g\n', sharpVig);
fprintf('  De-vignetted (AGBE)       : %.4g\n', sharpDevig);

figure('Position', [100 100 1400 600]);

subplot(2,3,1);
imagesc(origGray); axis image off;
colormap(gca, gray(256));
title('Original (IP, No Vignetting)', 'FontSize', 12, 'FontWeight', 'bold');
colorbar;

subplot(2,3,2);
imagesc(vigGray); axis image off;
colormap(gca, gray(256));
title('Vignetted (IP, Optical cos^4(\theta))', 'FontSize', 12, 'FontWeight', 'bold');
colorbar;

subplot(2,3,3);
imagesc(devigGray); axis image off;
colormap(gca, gray(256));
title(sprintf('De-Vignetted (AGBE)\nPSNR = %.2f dB', psnrVal), ...
      'FontSize', 12, 'FontWeight', 'bold');
colorbar;

subplot(2,3,4);
imagesc(vignetteMapEst); axis image off;
colormap(gca, parula(256));
title('Estimated Vignetting Map (AGBE)', 'FontSize', 12, 'FontWeight', 'bold');
colorbar;

subplot(2,3,5);
imagesc(residual); axis image off;
colormap(gca, jet(256));
title('Residual (De-Vig - Original)', 'FontSize', 12, 'FontWeight', 'bold');
colorbar;

subplot(2,3,6);
bar([sharpOrig, sharpVig, sharpDevig]);
set(gca, 'XTickLabel', {'Orig','Vig','De-Vig'});
ylabel('Var(Laplacian)');
title('Sharpness Comparison', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
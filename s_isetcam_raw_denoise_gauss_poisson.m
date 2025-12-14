ieInit;
ieSessionSet('waitbar', false);

img = imread('ss3.png');
scene = sceneFromFile(img, 'rgb', 100);

oi = oiCreate;
oi = oiCompute(oi, scene);

sensorBase = sensorCreate('bayer (rggb)');
sensorBase = sensorSet(sensorBase, 'fov', sceneGet(scene,'fov'), oi);
sensorBase = sensorSet(sensorBase, 'auto exposure', 'on');
sensorBase = sensorSet(sensorBase, 'noise flag', 0);

sensorClean = sensorCompute(sensorBase, oi);
rawClean = squeeze(sensorGet(sensorClean, 'volts'));

r = size(rawClean,1);
c = size(rawClean,2);
r = r - mod(r,2);
c = c - mod(c,2);
rawClean = rawClean(1:r,1:c);

rawScale  = max(rawClean(:));
rawCleanN = rawClean / max(rawScale, eps);

noiseType = 'gaussian';
noiseParams = struct();
noiseParams.sigma = 0.03;
noiseParams.photonsMax = 2000;

rawNoisyN = addControlledNoise(rawCleanN, noiseType, noiseParams);

sensorNoisy = sensorClean;
sensorNoisy = sensorSet(sensorNoisy, 'volts', rawNoisyN * rawScale);

denParamsG = struct('gaussSigma', 0.9);

sigmaEst01 = std(rawNoisyN(:) - rawCleanN(:));
denParamsB = struct('sigma', sigmaEst01);

rawMedN   = denoiseRawRGGB(rawNoisyN, 'median',   struct());
rawGaussN = denoiseRawRGGB(rawNoisyN, 'gaussian', denParamsG);
rawBm3dN  = denoiseRawRGGB(rawNoisyN, 'bm3d',     denParamsB);

sensorMed   = sensorNoisy; sensorMed   = sensorSet(sensorMed,   'volts', rawMedN   * rawScale);
sensorGauss = sensorNoisy; sensorGauss = sensorSet(sensorGauss, 'volts', rawGaussN * rawScale);
sensorBm3d  = sensorNoisy; sensorBm3d  = sensorSet(sensorBm3d,  'volts', rawBm3dN  * rawScale);

ip0 = ipCreate;
ip0 = ipSet(ip0, 'internal cs', 'XYZ');
ip0 = ipSet(ip0, 'conversion method sensor', 'MCC Optimized');
ip0 = ipSet(ip0, 'illuminant correction method', 'gray world');
ip0 = ipSet(ip0, 'demosaic method', 'adaptive laplacian');

ipClean = ipSet(ip0,'name','Clean');         ipClean = ipCompute(ipClean, sensorClean);
ipNoisy = ipSet(ip0,'name','Noisy');         ipNoisy = ipCompute(ipNoisy, sensorNoisy);
ipM     = ipSet(ip0,'name','Median pre');    ipM     = ipCompute(ipM,     sensorMed);
ipG     = ipSet(ip0,'name','Gaussian pre');  ipG     = ipCompute(ipG,     sensorGauss);
ipB     = ipSet(ip0,'name','BM3D pre');      ipB     = ipCompute(ipB,     sensorBm3d);

rgbClean = clip01(ipGet(ipClean,'srgb'));
rgbNoisy = clip01(ipGet(ipNoisy,'srgb'));
rgbM     = clip01(ipGet(ipM,    'srgb'));
rgbG     = clip01(ipGet(ipG,    'srgb'));
rgbB     = clip01(ipGet(ipB,    'srgb'));

mRawNoisy = metricPack(rawCleanN, rawNoisyN);
mRawMed   = metricPack(rawCleanN, rawMedN);
mRawG     = metricPack(rawCleanN, rawGaussN);
mRawB     = metricPack(rawCleanN, rawBm3dN);

mSrgbNoisy = metricPack(rgbClean, rgbNoisy);
mSrgbMed   = metricPack(rgbClean, rgbM);
mSrgbG     = metricPack(rgbClean, rgbG);
mSrgbB     = metricPack(rgbClean, rgbB);

fprintf('\n=== RAW metrics (vs clean RAW) ===\n');
printMetrics('RAW Noisy         ', mRawNoisy);
printMetrics('RAW Median        ', mRawMed);
printMetrics('RAW Gaussian      ', mRawG);
printMetrics('RAW BM3D    ', mRawB);

fprintf('\n=== sRGB metrics (vs clean sRGB) ===\n');
printMetrics('sRGB Noisy        ', mSrgbNoisy);
printMetrics('sRGB Median       ', mSrgbMed);
printMetrics('sRGB Gaussian     ', mSrgbG);
printMetrics('sRGB BM3D   ', mSrgbB);

g0 = wbGains(ipClean);
g1 = wbGains(ipNoisy);
g2 = wbGains(ipM);
g3 = wbGains(ipG);
g4 = wbGains(ipB);

fprintf('\n=== WB gains (normalized so G=1) [R G B] ===\n');
fprintf('Clean         : [%.4f  %.4f  %.4f]\n', g0(1), g0(2), g0(3));
fprintf('Noisy         : [%.4f  %.4f  %.4f]\n', g1(1), g1(2), g1(3));
fprintf('Median pre    : [%.4f  %.4f  %.4f]\n', g2(1), g2(2), g2(3));
fprintf('Gaussian pre  : [%.4f  %.4f  %.4f]\n', g3(1), g3(2), g3(3));
fprintf('BM3D pre      : [%.4f  %.4f  %.4f]\n', g4(1), g4(2), g4(3));

figure('Name','RAW mosaic (normalized)');
tiledlayout(2,5,'Padding','compact');
nexttile; imagesc(rawCleanN); axis image off; colormap gray; title('RAW clean');
nexttile; imagesc(rawNoisyN); axis image off; colormap gray; title('RAW noisy');
nexttile; imagesc(rawMedN);   axis image off; colormap gray; title('RAW median');
nexttile; imagesc(rawGaussN); axis image off; colormap gray; title('RAW gaussian');
nexttile; imagesc(rawBm3dN);  axis image off; colormap gray; title('RAW BM3D');
nexttile; imagesc(abs(rawNoisyN - rawCleanN)); axis image off; colormap gray; title('|noisy-clean|');
nexttile; imagesc(abs(rawMedN   - rawCleanN)); axis image off; colormap gray; title('|med-clean|');
nexttile; imagesc(abs(rawGaussN - rawCleanN)); axis image off; colormap gray; title('|gau-clean|');
nexttile; imagesc(abs(rawBm3dN  - rawCleanN)); axis image off; colormap gray; title('|bm3d-clean|');
nexttile; imagesc(abs(rawBm3dN  - rawMedN));   axis image off; colormap gray; title('|bm3d-med|');

figure('Name','sRGB after ISP (RAW denoise)');
tiledlayout(2,5,'Padding','compact');
nexttile; imshow(rgbClean); title('sRGB clean');
nexttile; imshow(rgbNoisy); title('sRGB noisy');
nexttile; imshow(rgbM);     title('sRGB median');
nexttile; imshow(rgbG);     title('sRGB gaussian');
nexttile; imshow(rgbB);     title('sRGB BM3D');
scaleDiff = 5;
nexttile; imshow(clip01(abs(rgbNoisy - rgbClean) * scaleDiff)); title('|noisy-clean| x5');
nexttile; imshow(clip01(abs(rgbM     - rgbClean) * scaleDiff)); title('|med-clean| x5');
nexttile; imshow(clip01(abs(rgbG     - rgbClean) * scaleDiff)); title('|gau-clean| x5');
nexttile; imshow(clip01(abs(rgbB     - rgbClean) * scaleDiff)); title('|bm3d-clean| x5');
nexttile; imshow(clip01(abs(rgbB     - rgbM    ) * scaleDiff)); title('|bm3d-med| x5');

function y = clip01(x)
y = min(max(double(x),0),1);
end

function noisy = addControlledNoise(clean, type, params)
sigma      = getOr(params, 'sigma',      0.02);
photonsMax = getOr(params, 'photonsMax', 2000);
clean = double(clean);
switch lower(type)
    case 'gaussian'
        noisy = clean + sigma * randn(size(clean));
        noisy = clip01(noisy);
    case 'poisson'
        lambda = clip01(clean) * photonsMax;
        if exist('poissrnd','file')
            counts = poissrnd(lambda);
        else
            counts = lambda + sqrt(lambda) .* randn(size(lambda));
        end
        noisy = counts / photonsMax;
        noisy = clip01(noisy);
    otherwise
        error('Unknown noise type');
end
end

function rawOut = denoiseRawRGGB(rawIn, method, params)
rawOut = rawIn;
R  = rawIn(1:2:end, 1:2:end);
G1 = rawIn(1:2:end, 2:2:end);
G2 = rawIn(2:2:end, 1:2:end);
B  = rawIn(2:2:end, 2:2:end);

switch lower(method)
    case 'bm3d'
        sigma01 = getOr(params, 'sigma', 0.01);
        try
            if ~exist('BM3D','file')
                error('BM3D_not_found');
            end
            sigma255 = sigma01 * 255;
            R  = bm3dGray01(R,  sigma255);
            G1 = bm3dGray01(G1, sigma255);
            G2 = bm3dGray01(G2, sigma255);
            B  = bm3dGray01(B,  sigma255);
        catch
            win = [5 5];
            noiseVar = sigma01^2;
            R  = wiener2(R,  win, noiseVar);
            G1 = wiener2(G1, win, noiseVar);
            G2 = wiener2(G2, win, noiseVar);
            B  = wiener2(B,  win, noiseVar);
        end
    otherwise
        f = makeDenoiseFcn(method, params);
        R  = f(R);  G1 = f(G1);  G2 = f(G2);  B  = f(B);
end

rawOut(1:2:end, 1:2:end) = R;
rawOut(1:2:end, 2:2:end) = G1;
rawOut(2:2:end, 1:2:end) = G2;
rawOut(2:2:end, 2:2:end) = B;
rawOut = clip01(rawOut);
end

function y = bm3dGray01(x01, sigma255)
z255 = clip01(x01) * 255;
[~, y255] = BM3D(1, z255, sigma255);
y = clip01(y255 / 255);
end

function f = makeDenoiseFcn(method, params)
switch lower(method)
    case 'gaussian'
        s = getOr(params, 'gaussSigma', 0.9);
        if exist('imgaussfilt','file')
            f = @(x) imgaussfilt(x, s);
        else
            h = fspecial('gaussian', max(3,2*ceil(3*s)+1), s);
            f = @(x) imfilter(x, h, 'replicate');
        end
    case 'median'
        f = @(x) medfilt2(x, [3 3]);
    otherwise
        error('Unknown denoise method');
end
end

function v = getOr(s, field, default)
if isstruct(s) && isfield(s, field)
    v = s.(field);
else
    v = default;
end
end

function m = metricPack(ref, test)
ref  = double(ref);
test = double(test);
d = test - ref;
m.MSE   = mean(d(:).^2);
m.MAE   = mean(abs(d(:)));
m.SNRdB = 20*log10(norm(ref(:)) / max(norm(d(:)), eps));
if exist('psnr','file')
    m.PSNR = psnr(test, ref);
else
    m.PSNR = NaN;
end
if ndims(ref)==3
    refY  = rgb2gray(clip01(ref));
    testY = rgb2gray(clip01(test));
else
    refY = ref; testY = test;
end
if exist('ssim','file')
    m.SSIM = ssim(testY, refY);
else
    m.SSIM = NaN;
end
if exist('niqe','file')
    m.NIQE = niqe(testY);
else
    m.NIQE = NaN;
end
if exist('corr2','file')
    m.Corr2 = corr2(testY, refY);
else
    m.Corr2 = NaN;
end
end

function printMetrics(label, m)
fprintf('%s | MSE: %.6g  MAE: %.6g  SNR(dB): %.2f  PSNR: %.2f  SSIM: %.4f  NIQE: %.3f  Corr2: %.4f\n', ...
    label, m.MSE, m.MAE, m.SNRdB, m.PSNR, m.SSIM, m.NIQE, m.Corr2);
end

function g = wbGains(ipObj)
T = [];
gotT = false;
try
    T = ipGet(ipObj, 'illuminant correction transform');
    gotT = true;
catch
end
if ~gotT
    try
        T = ipGet(ipObj, 'illuminant correction');
        gotT = true;
    catch
    end
end
if ~gotT
    g = [1; 1; 1];
    return;
end
g = ones(3,1);
if isnumeric(T) && all(size(T) == [3 3])
    g = diag(T);
end
if g(2) ~= 0
    g = g / g(2);
end
end

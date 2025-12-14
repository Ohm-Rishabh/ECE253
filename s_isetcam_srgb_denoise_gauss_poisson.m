%% s_isetcam_srgb_denoise_gauss_poisson.m
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
noiseParams.sigma      = 0.03;
noiseParams.photonsMax = 2000;

rawNoisyN = addControlledNoise(rawCleanN, noiseType, noiseParams);

sensorNoisy = sensorClean;
sensorNoisy = sensorSet(sensorNoisy, 'volts', rawNoisyN * rawScale);

ip0 = ipCreate;
ip0 = ipSet(ip0, 'internal cs', 'XYZ');
ip0 = ipSet(ip0, 'conversion method sensor', 'MCC Optimized');
ip0 = ipSet(ip0, 'illuminant correction method', 'gray world');
ip0 = ipSet(ip0, 'demosaic method', 'adaptive laplacian');

ipClean = ipSet(ip0,'name','Clean');  ipClean = ipCompute(ipClean, sensorClean);
ipNoisy = ipSet(ip0,'name','Noisy');  ipNoisy = ipCompute(ipNoisy, sensorNoisy);

rgbClean = clip01(ipGet(ipClean,'srgb'));
rgbNoisy = clip01(ipGet(ipNoisy,'srgb'));

denParamsG = struct('gaussSigma', 1.0);

sigmaEst01 = std(rgbNoisy(:) - rgbClean(:));
denParamsB = struct('sigma', sigmaEst01, 'wienerSize', [5 5]);

rgbMed = denoiseRgbImage(rgbNoisy, 'median',   struct());
rgbG   = denoiseRgbImage(rgbNoisy, 'gaussian', denParamsG);
rgbB   = denoiseRgbImage(rgbNoisy, 'bm3d',     denParamsB);

mRawNoisy = metricPack(rawCleanN, rawNoisyN);

fprintf('\n=== RAW metrics (noisy vs clean RAW) ===\n');
printMetrics('RAW Noisy        ', mRawNoisy);

mSrgbNoisy = metricPack(rgbClean, rgbNoisy);
mSrgbMed   = metricPack(rgbClean, rgbMed);
mSrgbG     = metricPack(rgbClean, rgbG);
mSrgbB     = metricPack(rgbClean, rgbB);

fprintf('\n=== sRGB metrics (vs clean sRGB) ===\n');
printMetrics('sRGB Noisy         ', mSrgbNoisy);
printMetrics('sRGB Median        ', mSrgbMed);
printMetrics('sRGB Gaussian      ', mSrgbG);
printMetrics('sRGB BM3D  ', mSrgbB);

g0 = wbGains(ipClean);
g1 = wbGains(ipNoisy);

fprintf('\n=== WB gains from ISP (normalized so G=1) [R G B] ===\n');
fprintf('Clean         : [%.4f  %.4f  %.4f]\n', g0(1), g0(2), g0(3));
fprintf('Noisy         : [%.4f  %.4f  %.4f]\n', g1(1), g1(2), g1(3));

figure('Name','RAW mosaic (normalized, controlled Gaussian/Poisson noise)');
tiledlayout(1,3,'Padding','compact');
nexttile; imagesc(rawCleanN); axis image off; colormap gray; title('RAW clean (original)');
nexttile; imagesc(rawNoisyN); axis image off; colormap gray; title('RAW noisy');
nexttile; imagesc(abs(rawNoisyN - rawCleanN)); axis image off; colormap gray; title('|noisy - clean|');

figure('Name','sRGB denoising AFTER ISP (late stage)');
tiledlayout(2,3,'Padding','compact');
nexttile; imshow(rgbClean); title('sRGB clean');
nexttile; imshow(rgbNoisy); title('sRGB noisy');
nexttile; imshow(rgbMed);   title('Median denoised');
nexttile; imshow(rgbG);     title('Gaussian denoised');
nexttile; imshow(rgbB);     title('BM3D');

function y = clip01(x)
y = min(max(double(x),0),1);
end

function noisy = addControlledNoise(clean, type, params)
if nargin < 3, params = struct(); end
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
        error('Unknown noise type: %s', type);
end
end

function rgbOut = denoiseRgbImage(rgbIn, method, params)
if nargin < 3, params = struct(); end
rgbIn = clip01(rgbIn);

switch lower(method)
    case 'bm3d'
        sigma01 = getOr(params, 'sigma', 0.01);
        win     = getOr(params, 'wienerSize', [5 5]);
        try
            if ~exist('CBM3D','file')
                error('CBM3D_not_found');
            end
            z255     = rgbIn * 255;
            sigma255 = sigma01 * 255;
            [~, y255] = CBM3D(1, z255, sigma255);
            rgbOut = clip01(y255 / 255);
        catch
            rgbOut = wienerRgb(rgbIn, win, sigma01^2);
        end

    otherwise
        rgbOut = rgbIn;
        for ch = 1:size(rgbIn,3)
            plane = rgbIn(:,:,ch);
            f = makeDenoiseFcn(method, params);
            rgbOut(:,:,ch) = f(plane);
        end
        rgbOut = clip01(rgbOut);
end
end

function out = wienerRgb(rgbIn, win, noiseVar)
out = rgbIn;
for ch = 1:size(rgbIn,3)
    out(:,:,ch) = wiener2(rgbIn(:,:,ch), win, noiseVar);
end
out = clip01(out);
end

function f = makeDenoiseFcn(method, params)
switch lower(method)
    case 'gaussian'
        s = getOr(params, 'gaussSigma', 1.0);
        if exist('imgaussfilt','file')
            f = @(x) imgaussfilt(x, s);
        else
            h = fspecial('gaussian', max(3,2*ceil(3*s)+1), s);
            f = @(x) imfilter(x, h, 'replicate');
        end
    case 'median'
        f = @(x) medfilt2(x, [3 3]);
    otherwise
        error('Unknown denoise method: %s', method);
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
    if exist('rgb2gray','file')
        refY  = rgb2gray(clip01(ref));
        testY = rgb2gray(clip01(test));
    else
        refY  = 0.2989*ref(:,:,1) + 0.5870*ref(:,:,2) + 0.1140*ref(:,:,3);
        testY = 0.2989*test(:,:,1) + 0.5870*test(:,:,2) + 0.1140*test(:,:,3);
        refY  = clip01(refY);
        testY = clip01(testY);
    end
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

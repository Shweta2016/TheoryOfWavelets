%Read image
bits = 8;
cutoff = 0.85;
%cutoff = 0.90;
%cutoff = 0.95;
len = 256;
A = imread('image.jpg');
A = rgb2gray(A);
A = double(A);
A = imresize(A,[len,len],'bicubic');
figure(1);
imagesc(A);
colormap gray(256);
title('Original Image');

%D6 Encoding
[bins,sign,c] = d6_enc_log(A,bits,cutoff);
original_bytes = len^2;
% Apply gunzip
ratio = compress_lossless(bins,sign,original_bytes);
clear bins sign
[bins,sign] = uncompress_lossless('Bins6','Signs6');
sign(sign==0) = -1;
% Dequantize and reconstruct image
A_rc = inv_d4_dequant(bins,sign,c);
figure(3);
imagesc(A_rc);
colormap gray(256);
title('Reconstructed Image');

%PSNR
psnr = peak_signal_noise_ratio(A,A_rc);
fprintf('Compression ratio: %f\n', ratio);
fprintf('Peak signal to noise ratio cutoff 0.85: %f\n', psnr);

function mse = mean_squared_error(I,K)
    [m,n] = size(I);
    mse = 0;
    for i = 1:m
    for j = 1:n
    mse = mse + ( I(i,j) - K(i,j) )^2;
    end
    end
    mse = mse / (m*n);
end

function psnr = peak_signal_noise_ratio(I,K)
    mse = mean_squared_error(I,K);
    psnr = 10*log10( ( max(I(:))^2 )/mse );
end

function ratio = compress_lossless(bins,sgn,original_bytes)
    working_path = pwd;
    % write array bins to file Bins
    FILE1='Bins6';
    fid=fopen(FILE1,'w');
    count=fwrite(fid,bins);
    status=fclose(fid);
    % write array sgn to file Sign
    FILE2='Signs6';
    fid=fopen(FILE2,'w');
    count=fwrite(fid,sgn);
    status=fclose(fid);
    % apply gzip
    gzip(FILE1);
    gzip(FILE2);
    % Number of bytes of Bins after gzip.
    FILE1_BYTES=strcat(working_path,'/',FILE1,'.gz');
    s=dir(FILE1_BYTES);
    compressed1_bytes = s.bytes;
    % Number of bytes of Sign after gzip.
    FILE2_BYTES=strcat(working_path,'/',FILE2,'.gz');
    s=dir(FILE2_BYTES);
    compressed2_bytes = s.bytes;
    % Compression ratio
    ratio = original_bytes/(compressed1_bytes+compressed2_bytes);
end

function [bins,sgn] = uncompress_lossless(FILE1,FILE2)
    working_path = pwd;
    GZIP1=strcat(FILE1,'.gz');
    GZIP2=strcat(FILE2,'.gz');
    gunzip(GZIP1);
    gunzip(GZIP2);
    fid=fopen(FILE1,'r','l');
    bins=fread(fid);
    status=fclose(fid);
    fid=fopen(FILE2,'r','l');
    sgn=fread(fid);
    status=fclose(fid);
end

% Encoding Function
function [y,s,c] = d6_enc_log(B,bits,cutoff)
    len = length(B);
    % build Haar filter matrix
    Q=[1 1;1 -1];
    I = eye(len);
    H = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
    % D4 Encoding Matrix
    %Find coefficients
    syms h0 h1 h2 h3 h4 h5
    S=solve(h0^2+h1^2+h2^2+h3^2+h4^2+h5^2==1,...
        h0+h1+h2+h3+h4+h5==sqrt(2),...
        h0-h1+h2-h3+h4-h5==0,...
        h0*h2+h1*h3+h2*h4+h3*h5==0,...
        h0*h4+h1*h5==0,...
        h1-2*h2+3*h3-4*h4+5*h5==0,...
        -h1+4*h2-9*h3+16*h4-25*h5==0);


    h0=double(S.h0);
    h1=double(S.h1);
    h2=double(S.h2);
    h3=double(S.h3);
    h4=double(S.h4);
    h5=double(S.h5);
    %Symbol hk
    h = zeros(1,6);
    h(1) = h0(1);
    h(2) = h1(1);
    h(3) = h2(1);
    h(4) = h3(1);
    h(5) = h4(1);
    h(6) = h5(1);
    Q1 = double([h(1) h(2); h(6) -h(5)]);
    Q2 = double([h(3) h(4); h(4) -h(3)]);
    Q3 = double([h(5) h(6); h(2) -h(1)]);
    I = eye(len);
    % build permutation matrix
    PT = I([1:2:len],:);
    PB = I([2:2:len],:);
    % encode image
    len = 256;
    for j = 1:log2(len)-2
        P = [PT(1:len/2,1:len);
        PB(1:len/2,1:len)];
        d1 = kron(I(1:len/2,1:len/2),Q1);
        d2 = kron(I(1:len/2,1:len/2),Q2);
        d3 = kron(I(1:len/2,1:len/2),Q3);
        D = d1 + circshift(d2,2,2) + circshift(d3,4,2);
        B(1:len,1:len)=P*D*B(1:len,1:len)*D'*P';
        len = len/2;
    end
    
    figure(2);
    image(B);
    colormap gray(16);
    title('Encoded Image');
    
    % apply thresholding
    len = 256;
    [len,len] = size(B);
    X = sort(abs(B(:)));
    th = X(floor(cutoff*len^2));
    BT = B;
    BT(abs(BT) < th) = 0;
    % Log Quantization
    x = BT(:);
    NP = 2^bits;
    NX = size(x);
    k = 1;
    a = abs(x);
    for n=1:NX
    if a(n) > th-eps;
    s(k) = sign(x(n));
    k=k+1;
    end
    end
    MX = max(a);
    c = zeros(NP,1);
    p = zeros(NP-1,1);
    c(1) = 0.;
    c(NP) = MX;
    p(1) = th;
    d = (MX/th)^(1/(NP-1));
    for n=2:NP-1
    p(n) = th*d^n;
    c(n) = (p(n-1)+p(n))/2;
    end
    % p(NP-1) = p(NP-1)-eps;
    p(NP-1) = p(NP-1)+eps;
    y = quantiz(a,p);
end

% Decoding Function
function A = inv_d6_dequant(y,s,c)
    BQ = c(y(:)+1);
    k=1;
    len = length(c);
    for n = 1:len^2
    if BQ(n) ~= 0
    BQ(n) = BQ(n)*s(k);
    k = k+1;
    end
    end
    BQ = reshape(BQ,[len,len]);
    len = length(BQ);
    % D4 Encoding Matrix
    %Find coefficients
    syms h0 h1 h2 h3 h4 h5
    S=solve(h0^2+h1^2+h2^2+h3^2+h4^2+h5^2==1,...
        h0+h1+h2+h3+h4+h5==sqrt(2),...
        h0-h1+h2-h3+h4-h5==0,...
        h0*h2+h1*h3+h2*h4+h3*h5==0,...
        h0*h4+h1*h5==0,...
        h1-2*h2+3*h3-4*h4+5*h5==0,...
        -h1+4*h2-9*h3+16*h4-25*h5==0);

    h0=double(S.h0);
    h1=double(S.h1);
    h2=double(S.h2);
    h3=double(S.h3);
    h4=double(S.h4);
    h5=double(S.h5);
    %Symbol hk
    
    h = zeros(1,6);
    h(1) = h0(1);
    h(2) = h1(1);
    h(3) = h2(1);
    h(4) = h3(1);
    h(5) = h4(1);
    h(6) = h5(1);
    Q1 = double([h(1) h(2); h(6) -h(5)]);
    Q2 = double([h(3) h(4); h(4) -h(3)]);
    Q3 = double([h(5) h(6); h(2) -h(1)]);
    I = eye(len);
    A = BQ;

    for len = 2.^(2:log2(len))
    % build permutation matrix
    PT = I([1:2:len],:);
    PB = I([2:2:len],:);
    P = [PT(1:len/2,1:len);PB(1:len/2,1:len)];
    d1 = kron(I(1:len/2,1:len/2),Q1);
    d2 = kron(I(1:len/2,1:len/2),Q2);
    d3 = kron(I(1:len/2,1:len/2),Q3);
    D = d1 + circshift(d2,2,2) + circshift(d3,4,2);
    A(1:len,1:len) = D'*P'*A(1:len,1:len)*P*D;
    end
end

% Decoding Function
function A = inv_d4_dequant(y,s,c)
    BQ = c(y(:)+1);
    k=1;
    len = length(c);
    for n = 1:len^2
    if BQ(n) ~= 0
    BQ(n) = BQ(n)*s(k);
    k = k+1;
    end
    end
    BQ = reshape(BQ,[len,len]);
    len = length(BQ);
    % D4 Encoding Matrix
    syms h0 h1 h2 h3 h4 h5
    S=solve(h0^2+h1^2+h2^2+h3^2+h4^2+h5^2==1,...
        h0+h1+h2+h3+h4+h5==sqrt(2),...
        h0-h1+h2-h3+h4-h5==0,...
        h0*h2+h1*h3+h2*h4+h3*h5==0,...
        h0*h4+h1*h5==0,...
        h1-2*h2+3*h3-4*h4+5*h5==0,...
        -h1+4*h2-9*h3+16*h4-25*h5==0);

    h0=double(S.h0);
    h1=double(S.h1);
    h2=double(S.h2);
    h3=double(S.h3);
    h4=double(S.h4);
    h5=double(S.h5);
    %Symbol hk
    
    h = zeros(1,6);
    h(1) = h0(1);
    h(2) = h1(1);
    h(3) = h2(1);
    h(4) = h3(1);
    h(5) = h4(1);
    h(6) = h5(1);
    Q1 = double([h(1) h(2); h(6) -h(5)]);
    Q2 = double([h(3) h(4); h(4) -h(3)]);
    Q3 = double([h(5) h(6); h(2) -h(1)]);
    I = eye(len);
    A = BQ;
    % Decode image
    for len = 2.^(3:log2(len))
    % build permutation matrix
    PT = I([1:2:len],:);
    PB = I([2:2:len],:);
    P = [PT(1:len/2,1:len);
    PB(1:len/2,1:len)];
    d1 = kron(I(1:len/2,1:len/2),Q1);
    d2 = kron(I(1:len/2,1:len/2),Q2);
    d3 = kron(I(1:len/2,1:len/2),Q3);
    D = d1 + circshift(d2,2,2) + circshift(d3,4,2);
    A(1:len,1:len) = D'*P'*A(1:len,1:len)*P*D;
    end
end


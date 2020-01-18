%read image
A=imread('image.jpg');A=rgb2gray(A);
A=double(A);
len=256;B=imresize(A,[len,len],'bicubic');
figure(1);imagesc(B);colormap gray(256);
title('Original Image')
%build Haar filter matrix
Q=[1 1;1 -1];I = eye(len);
H = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
%build permutation matrix
PT = I([1:2:len],:);PB = I([2:2:len],:);
%encode image
for j = 1:log2(len)
    P = [PT(1:len/2,1:len); PB(1:len/2,1:len)];
    H = H(1:len,1:len);
    B(1:len,1:len)=P*H*B(1:len,1:len)*H'*P';
    len = len/2;
end
figure(2);image(B);colormap gray(16);BB=B;
title('Haar Encoded Image')
% Get threshold
cutoff = 0.85;  % Output with cutoff=.85,.90 and .95
[len,len] = size(B);
X = sort(abs(B(:)));
th = X(floor(cutoff*len^2));BQ=B;BQ(abs(B)<=th)=0;
% Thresholding and log quantization
% Inputs
% x: 2-D array
x = B(:); % Tried with image A
% th: threshold
% bits: number of bits for quantization
bits = 8;
% Outputs:
% y: ?bins? (thresholded and quantized abs(x))
% s: sign of nonzero x post thresholding
% c: codebook
NP = 2^bits;
NX=length(x);
k=1;
a=abs(x(:));
for n=1:NX
    if a(n)> th; 
        s(k)=sign(x(n));
        k=k+1;
    end
end
s1 = s;
s=s';
MX = max(a);
c=zeros(NP,1);
p=zeros(NP-1,1);
c(1)=0.;
c(NP)=MX;
p(1)=th;
d =(MX/th)^(1/(NP-1));
for n=2:NP-1
    p(n)=th*d^n; 
    c(n)= (p(n-1)+p(n))/2;
end
p(NP-1)=p(NP-1)-eps;
y = quantiz(a,p);
y1 = reshape(y,[256,256]);
figure(3);image(reshape(y,[256,256]));colormap gray(16);
title('Quantized with Cutoff 85%')
%Dequantization
BQ = c(y(:)+1);
k=1;
len = length(y);
s=s(:);BQ=BQ(:);
for n=1:65536
    if BQ(n)~= 0
        BQ(n)=BQ(n)*s(k);
        k=k+1;
    end
end
BQ=reshape(BQ,[256,256]);
figure(4);image(BQ);colormap gray(16);BBB=BQ;
title('DeQuantized with Cutoff 85%')
%Decode image
len = 256;
len1 = 2;
Q1=[1 1;1 -1];I1 = eye(len);
H1 = kron(I1(1:len/2,1:len/2),Q1)/sqrt(2);
for j = 1:log2(len)
    P1 = [PT(1:len1/2,1:len1); PB(1:len1/2,1:len1)];
    Q1=[1 1;1 -1];I1 = eye(len1);
    H1 = kron(I1(1:len1/2,1:len1/2),Q1)/sqrt(2);
    BQ(1:len1,1:len1)=H1'*P1'*BQ(1:len1,1:len1)*P1*H1;
    len1 = len1*2;
end
figure(5);image(BQ);colormap gray(256);
title('Decoded with Cutoff 85%')





working_path = pwd;


FILE1='Bins';fid=fopen(FILE1,'w');count=fwrite(fid,y1);status=fclose(fid);
%write array sgn to file Sign
FILE2='Sign';fid=fopen(FILE2,'w');count=fwrite(fid,s1);status=fclose(fid);
%apply gzip
gzip(FILE1);gzip(FILE2);
original_bytes = 65536+9831;

%Number of bytes of Bins after gzip.
FILE1_BYTES=strcat(working_path,'/',FILE1,'.gz');
s=dir(FILE1_BYTES);compressed1_bytes = s.bytes;
%Number of bytes of Sign after gzip.
FILE2_BYTES=strcat(working_path,'/',FILE2,'.gz');
s=dir(FILE2_BYTES);compressed2_bytes = s.bytes;
%Compression ratio
ratio = original_bytes/(compressed1_bytes+compressed2_bytes);
display(original_bytes);
display(compressed1_bytes);
display(compressed2_bytes);
display(ratio);

working_path = pwd;
GZIP1=strcat('Bins','.gz');GZIP2=strcat('Sign','.gz');
gunzip(GZIP1);gunzip(GZIP2);
fid=fopen('Bins','r','l');bins=fread(fid);status=fclose(fid);
fid=fopen('Sign','r','l');sgn=fread(fid);status=fclose(fid);


%Peak to noise ratio
% Original Image = A and Reconstructed Image = BQ
MSE = 0;
MSE1 = 0;
for i = 1:255
    for j = 1:255
        val = A(i,j) - y1(i,j);
        sq = val * val;
        MSE1 = MSE1 + sq;
    end
end
mn = 256*256;
MSE = MSE1 / mn;

display(MSE);

%PSNR : Peak Signal to Noise Ratio
PSNR = 0;
PSNR = (20 * log10(255)) - (10 * log10(MSE));
display(PSNR);

        
        
        
        
        
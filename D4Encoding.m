%read image
A=imread('image.jpg');A=rgb2gray(A);
A=double(A);
len=256;B=imresize(A,[len,len],'bicubic');
figure(1);imagesc(B);colormap gray(256);
title('Original Image');daspect([1 1 1]);
Original = B;
I = eye(len);
%build D4 filter matrix
h0 = 0.4830;
h1 = 0.8365;
h2 = 0.2241;
h3 = -0.1294;
Q=[h0 h1 h2 h3;h3 -h2 h1 -h0];
H1 = zeros(2,4);
H1(1:2,1:4) = Q;
r=1;
c=1;
H = zeros(2,4);
for ll = 1:127
    H(c:c+1,r:r+3) = H1;
    c=c+2;
    r=r+2;
end
H(255:256,255:256) = [h0 h1;h3 -h2];
H(255:256,1:2) = [h2 h3;h1 -h0];
H2 = H;
%disp(H);
%size1 = size(H);
%disp(size1);

%build permutation matrix
PT = I([1:2:len],:);PB = I([2:2:len],:);
%encode image
for j = 1:7
    P = [PT(1:len/2,1:len); PB(1:len/2,1:len)];
    H = H(1:len,1:len);
    %disp(size(P));
    %disp(size(H));
    %disp(size(B));
    B(1:len,1:len)=P*H*B(1:len,1:len)*H'*P';
    len = len/2;
end
figure(2);image(B);colormap gray(16);BB=B;
title('D4 Encoded Image')
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
title('Quantized with Cutoff %')
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
title('DeQuantized with Cutoff %');
%Decode image
len = 256;
len1 = 2;
len2 = 4;
for j = 1:7
    P1 = [PT(1:len2/2,1:len2); PB(1:len2/2,1:len2)];
    H3 = H2(1:len2,1:len2);
    %disp(size(P1));
    %disp(size(H3));
    %disp(size(BQ));
    BQ(1:len2,1:len2)=H3'*P1'*BQ(1:len2,1:len2)*P1*H3;
    len2 = len2*2;
    
    
end
figure(5);image(BQ);colormap gray(256);
title('Decoded with Cutoff %');daspect([1 1 1]);
working_path = pwd;


FILE1='Bins1';fid=fopen(FILE1,'w');count=fwrite(fid,y1);status=fclose(fid);
%write array sgn to file Sign
FILE2='Sign1';fid=fopen(FILE2,'w');count=fwrite(fid,s1);status=fclose(fid);
%apply gzip
gzip(FILE1);gzip(FILE2);
original_bytes = 256^2;

%Number of bytes of Bins after gzip.
FILE1_BYTES=strcat(working_path,'/',FILE1,'.gz');
s=dir(FILE1_BYTES);compressed1_bytes = s.bytes;
%Number of bytes of Sign after gzip.
FILE2_BYTES=strcat(working_path,'/',FILE2,'.gz');
s=dir(FILE2_BYTES);compressed2_bytes = s.bytes;
%Compression ratio
ratio = original_bytes/(compressed1_bytes+compressed2_bytes);


working_path = pwd;
GZIP1=strcat('Bins1','.gz');GZIP2=strcat('Sign1','.gz');
gunzip(GZIP1);gunzip(GZIP2);
fid=fopen('Bins1','r','l');bins=fread(fid);status=fclose(fid);
fid=fopen('Sign1','r','l');sgn=fread(fid);status=fclose(fid);


%Peak to noise ratio
% Original Image = A and Reconstructed Image = BQ
npts=256^2; mse=sum((Original(:)-BQ(:)).^2)/npts; psnr=10*log10(255^2/mse);
fprintf('Haar Comp Ratio = %.4f \n',ratio)
fprintf('Haar PSNR = %.4f dB \n',psnr)

        
        
        
        
        
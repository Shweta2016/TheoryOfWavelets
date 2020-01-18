%Image Construction
X = ones(256,256,3);

X(25:225,25:26,:) = 0;
X(25:225,224:225,:) = 0;
X(25:26,25:225,:) = 0;
X(224:225,25:225,:) = 0;

for c = 25:224
    for r = 25:224
        if r == c
            X(r:r+1,c:c+1,:) = 0;
        end
        if (r+c == 249)
            X(r:r+1,c:c+1,:) = 0;
        end
          
    end
end
figure(1);image(X);
X=rgb2gray(X);
X = double(X);
len=256;
X=imresize(X,[len,len],'bicubic');
figure(2);imagesc(X);colormap gray(256);
title('Original Image')
%build Haar filter matrix
Q=[1 1;1 -1];I = eye(len);
H = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
%build permutation matrix
PT = I([1:2:len],:);PB = I([2:2:len],:);
P = [PT;PB];
B = P*H*X*H'*P';
figure(3);imagesc(B);colormap gray(16);
title('Haar Level 1 Encoded Image')

%encode image
for j = 1:log2(len)
    P = [PT(1:len/2,1:len); PB(1:len/2,1:len)];
    H = H(1:len,1:len);
    BB(1:len,1:len)=P*H*X(1:len,1:len)*H'*P';
    len = len/2;
end
figure(4);imagesc(BB);colormap gray(16);
title('Haar Encoded Image')





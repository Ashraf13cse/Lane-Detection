
function compress = dwt_db1(X,r,c)
% X=imread('f00001.png');
X=X(1:r,1:c);

[a1,a2]=size(X);

[cA1,cH1,cV1,cD1] = dwt2(X,'db1');
[C,S] = wavedec2(X,1,'db1');
A1 = wrcoef2('a',C,S,'db1',1);
H1 = wrcoef2('h',C,S,'db1',1); 
V1 = wrcoef2('v',C,S,'db1',1); 
D1 = wrcoef2('d',C,S,'db1',1); 


re_ima1 = idwt2(cA1,cH1,cV1,cD1,'db1'); 
re_ima=uint8(re_ima1);

% [C,S] = wavedec2(X,2,'db1');
% A2 = wrcoef2('a',C,S,'db1',2);
% A1 = wrcoef2('a',C,S,'db1',1);
% H1 = wrcoef2('h',C,S,'db1',1); 
% V1 = wrcoef2('v',C,S,'db1',1); 
% D1 = wrcoef2('d',C,S,'db1',1); 
% H2 = wrcoef2('h',C,S,'db1',2);
% V2 = wrcoef2('v',C,S,'db1',2); 
% D2 = wrcoef2('d',C,S,'db1',2);

% dec2d = [A2,A1,H1,V1,D1,H2,V2,D2];
dec2d = [A1,H1,V1,D1];
re_ima1 = waverec2(C,S,'db1'); 
re_ima=uint8(re_ima1);
n=1;

X=X(1:r,1:c);
X=double(X)-128;
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('db1');
[c,s]=wavedec2(uint8(X),n,Lo_D,Hi_D);

[thr,nkeep] = wdcbm2(uint8(dec2d),1.5,prod(s(1,:)));
[THR,SORH,KEEPAPP,CRIT] = ddencmp('cmp','wp',uint8(X));
[XC,TREED,PERF0,PERFL2] =wpdencmp(X,SORH,2,'db1',CRIT,THR,KEEPAPP);

XC=double(X)+128;
XO=uint8(XC);

compress=XO;
imshow(compress)
end

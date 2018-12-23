
clc
close all
clear all
for imageno=1:1
    if imageno<=9
A=imread('f00001.png');
%          A= imread(['f0000',num2str(imageno),'.png']);
%     elseif imageno<=99
%         A= imread(['f000',num2str(imageno),'.png']);
%     elseif imageno<=999
%         A= imread(['f00',num2str(imageno),'.png']);
    end
     figure,imshow(A),figure

B=A;    
[q w]=size(A);
w=w/3;
A=A(1:320,50:540,:);
fll=1;
ID=mean(A,2);
DD=min(ID(:,:,:));
[ro co]= find(ID == DD(:,:,3));

A=A(ro+40:end,1:end,:);
figure,imshow(A),figure
R=A;
%%

[r c s]=size(A);
AA=dwt_db(R,r,c);

    Z=zeros(r,c);
    for i=1:r
        for j=1:c
            Z(i,j)=max([A(i,j,1) A(i,j,2) A(i,j,3)]);
        end
    end
    
    
    kk=8;
    D=zeros(r,c);
    
    C=double(A);
    for i=1:r
        for j=kk+1:c-kk
            D(i,j)=(C(i,j)^2-C(i,j+kk)*C(i,j-kk));
        end
    end
    F=zeros(r,c);
   figure,imshow(uint32(D),[min(min(D)) max(max(D))])
    thval=prctile(reshape(D,1,size(D,1)*size(D,2)),95);
    for i=1:r
        for j=1:c
            if D(i,j)>thval
                F(i,j)=1;
            else
                F(i,j)=0;
            end
        end
    end
    
     figure,imshow(F)
    ZZ=ones(r,c);
fct=2;
    for iiii=1:3*round(r/4)
        ZZ(3*round(r/4)-iiii+1,:)=[zeros(1,round(fct*iiii)) ones(1,c-2*round(fct*iiii)) zeros(1,round(fct*iiii))];
    end
%     fact=71;
%     fact2=120;
%     for iiii=1:r
%         ZZ(r-iiii+1,:)=[zeros(1,fact2) ones(1,c-(fact+fact2)) zeros(1,fact)];
%         fact=fact+1;
%         fact2 = fact2+1;
%     end
% figure,imshow(ZZ),figure
    F=ZZ.*F;
%     F=bwareaopen(F,10);
    F=bwareaopen(F,30);
    [labeledImage, numberOfRegions] = bwlabel(F,8);
      figure,imshow(labeledImage)
%    labeledImage = bwareaopen(labeledImage, 20);
%     imshow(labeledImage);
    fprintf('Number of regions = %d\n', numberOfRegions)
    coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
     figure,imshow(coloredLabels);
%     labeledImage=bwmorph(labeledImage,'skel');
    
    %
    measurements = regionprops(labeledImage, 'MajorAxisLength', 'PixelList')
    
    singleColumnIndexes = [];
    poss=[];

 xy=[];
 
    for k = 1 : numberOfRegions
       
        xy =  measurements(k).PixelList;
       
        x = xy(:, 1);
        a(k)=x(1);
        
        y = xy(:, 2);
        b(k)=y(1);
        length(y);
        

            
            CC = polyfit(x,y,1);
            %dal=CC(1);
            slope(k)=CC(1);
            constant(k)=round(CC(2)*10)/10
%             constant(k)
            
            
            
            Y = polyval(CC, 1);
%             YY=20;

      
    end
    slope=round(slope*100)/100;
%     slope=sort(slope,'descend');
%     constant=sort(constant,'descend');
%% start Clustering
cnt=1;
equalSlope=[];
vis=[];
for i=1:numberOfRegions
    vis(i)=0;
end
    for i=1:numberOfRegions
        m=slope(i);
        if(vis(i)==0)
             equalSlope(i)=cnt;
            for j=i+1:numberOfRegions
                if(abs(slope(i)- slope(j))<=.9 && vis(j)==0)
                   vis(j)=1;
                    equalSlope(j)=cnt;
                end      
            end
            cnt=cnt+1;
        end
        
    end
    
 
    

%% determine the cluster length
clusterMaxValue=[];
for i=1:numberOfRegions
    vis(i)=0;
    clusterMaxValue(equalSlope(i))=-1;
end

cnt=1;
sameLineSize=[];
Length=[];
sameLineIdx=[];
for ii=1:numberOfRegions
sameLineIdx(ii)=0;
end

for i=1:numberOfRegions
     m = equalSlope(i);
    % sameLineIdx(i)=i;
         xy =  measurements(i).PixelList;
         xx = xy(:,1);
         Length(i)=length(xx);
         M = slope(i);
         C = constant(i);
         %fprintf('Same Line\n');
        sz=0;
    Ln=[];
    for ii=1:numberOfRegions
        Ln(ii)=0;
    end
      for j=1:numberOfRegions
            if(equalSlope(i)==equalSlope(j) && i~=j)
               
                xy =  measurements(j).PixelList;
                x = xy(:, 1);
                y = xy(:, 2);
               flag=0;
                for k=1:length(x) 
                  if(abs(y(k) - (M*x(k) + C))<=1.0)
                    flag=1;
                    break;
                  end
                end 
                if flag==1
                    Ln(j)= 1;
                    sz=sz+length(x);
                end
                %fprintf('-------------\n');
                
            end    
      end
      if(clusterMaxValue(equalSlope(i)) < length(xx) + sz)
          
          clusterMaxValue(equalSlope(i)) = length(xx) + sz;
          sameLineIdx(i)=i;
          for k=1:numberOfRegions
              if Ln(k)==1
                  sameLineIdx(k)=i;
              end
          end
          
      end
     sameLineSize(i) = length(xx) + sz;
end 
%% end the determination of cluster size %%

temp= sort(sameLineSize,'descend');



greaterThesholdValueIdx=[];
c=1;
for i=1:numberOfRegions
    if sameLineSize(i)>20
        greaterThesholdValueIdx(c) = i;
        c=c+1;
    end
end
for i=1:numberOfRegions
    vis(i)=0;
end
fprintf("Count........... c");
% c


for i=1:c-1
    t=greaterThesholdValueIdx(i);
   for j=1:numberOfRegions
       if t == sameLineIdx(j)
           vis(j)=1;
       end
   end
end

poss=[];
 NegSlopMaxX=-1;
  NegSlopMaxY=-1;
   NegSlopMinX=10007;
    NegSlopMinY=10007;
    
     PosSlopMaxX=-1;
  PosSlopMaxY=-1;
   PosSlopMinX=10007;
    PosSlopMinY=10007;
    
for i=1:numberOfRegions
           xy =  measurements(i).PixelList;
            x = xy(:,1);
            y = xy(:,2);

            if vis(i)==0
                for p=1:length(y)
                    labeledImage(y(p),x(p))=0;
                end
            end
end
        
 for i=1:numberOfRegions
      xy =  measurements(i).PixelList;
      x = xy(:,1);
      y = xy(:,2);

     if vis(i)==1 && length(x) > 150
%          if(numberOfRegions<4)
            
            CC = polyfit(x,y,1)
            Y = polyval(CC, 1);
            YY=10;
            XX=round((YY-CC(2))/CC(1));
            YY2=300;
            XX2=round((YY2-CC(2))/CC(1));
            %% try 0 and 1 for full-line or fitted line
%             fullline=0;
            fullline=1;
            %%
            if fullline==1
                poss=[poss;XX YY XX2 YY2]  ;% for line start to end
            else
                if CC(1)<0
                    poss=[poss;min(x) max(y) max(x) min(y)];
                    
                else
                    poss=[poss;min(x) min(y) max(x) max(y)];
                    
                end
            end
%          end
     end
%             x1 = x(1);
%             y1 = y(1);
%             x2 = x(length(x));
%             y2 = y(length(x));
%             m = (y1-y2)/(x1-x2);
%            if m<0
%                     %poss=[poss;min(x) max(y) max(x) min(y)];
%                     NegSlopMaxX = max(NegSlopMaxX,max(x));
%                     NegSlopMaxY = max(NegSlopMaxY,max(y));
%                     NegSlopMinX = min(NegSlopMinX,min(x));
%                     NegSlopMinY = min(NegSlopMinY,min(y));
%            
%             else
%                     %poss=[poss;min(x) min(y) max(x) max(y)];
%                     PosSlopMaxX = max(PosSlopMaxX,max(x));
%                     PosSlopMaxY = max(PosSlopMaxY,max(y));
%                     PosSlopMinX = min(PosSlopMinX,min(x));
%                     PosSlopMinY = min(PosSlopMinY,min(y));
%            end
%             
%      end
     
     end
% if NegSlopMinX ~=10007 && NegSlopMaxY ~= -1 && NegSlopMaxX~= -1 && NegSlopMinY ~=10007
%     poss=[poss;NegSlopMinX NegSlopMaxY NegSlopMaxX NegSlopMinY];
% end
% if PosSlopMinX ~=10007 && PosSlopMinY ~= 10007 && PosSlopMaxX~= -1 && PosSlopMaxY~= -1
%  poss=[poss;PosSlopMinX PosSlopMinY PosSlopMaxX PosSlopMaxY];
% end


%sameLine = sort(sameLine,'descend');
     
    [labeledImage, numberOfRegions] = bwlabel(labeledImage,8);
    figure,imshow(labeledImage)
    fprintf('Number of lanes = %d\n', numberOfRegions);
    coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels

     II=insertShape(R,'line',poss,'Color','blue','LineWidth',5);
     figure,imshow(II)
%      II=insertShape(B,'line',poss2,'Color','blue','LineWidth',5);
%      figure,imshow(II)
toc
end


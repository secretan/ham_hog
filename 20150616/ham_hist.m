%% ham直方图统计

%function descriptor = ham_hist(InImage, xstart, ystart, xoffset, yoffset)
function descriptors_ = ham_hist(InImage, cWinx, cWiny, cBlkx, cBlky, cCelx, cCely, cSrdx, cSrdy);
descriptors = zeros((round((cWinx-cBlkx)/cSrdx)+1)*(round((cWiny-cBlky)/cSrdy)+1)*(cBlkx/cCellx)*(cBlky/cCelly),9);
[m,n] = size(InImage);
ZImage = zeros(m+2,n+2);
for i=2:m+1
    for j=2:n+1
        ZImage(i,j) = InImage(i-1,j-1);
    end
end
% 第一行==第二行
ZImage(1,:) = ZImage(2,:);
% 倒数第一行==倒数第二行
ZImage(m+2,:) = ZImage(m+1,:);
% 第一列==第二列
ZImage(:,1) = ZImage(:,2);
% 倒数第一列==倒数第二列
ZImage(:,n+2) = ZImage(:,n+1);

Hist = zeros(1,9);
for i = 1:xoffset
    for j=1:yoffset
        RectImage(i,j) = InImage(xstart+i,ystart+j);
    end
end

for i = 1:xoffset
    for j = 1:yoffset
        if RectImage(i,j) >= 0 && RectImage(i,j) < 20
            Hist(1,1) = Hist(1,1) + 1;
        elseif RectImage(i,j) >= 20 && RectImage(i,j) < 40
            Hist(1,2) = Hist(1,2) + 1;
        elseif RectImage(i,j) >= 40 && RectImage(i,j) < 60
            Hist(1,3) = Hist(1,3) + 1;            
        elseif RectImage(i,j) >= 60 && RectImage(i,j) < 80
            Hist(1,4) = Hist(1,4) + 1;
        elseif RectImage(i,j) >= 80 && RectImage(i,j) < 90
            Hist(1,5) = Hist(1,5) + 1;
        elseif RectImage(i,j) >= -90 && RectImage(i,j) < -80
            Hist(1,5) = Hist(1,5) + 1;
        elseif RectImage(i,j) >= -80 && RectImage(i,j) < -60
            Hist(1,6) = Hist(1,6) + 1;
        elseif RectImage(i,j) >= -60 && RectImage(i,j) < -40
            Hist(1,7) = Hist(1,7) + 1;
        elseif RectImage(i,j) >= -40 && RectImage(i,j) < -20      
            Hist(1,8) = Hist(1,8) + 1;      
        elseif RectImage(i,j) >= -20 && RectImage(i,j) < 0
            Hist(1,9) = Hist(1,9) + 1;
        end
    end
end
descriptor = Hist;

end
%% ham直方图统计

%function descriptor = ham_hist(InImage, xstart, ystart, xoffset, yoffset)
function descriptors_ = ham_hist2(Grad,Angle, cWinx, cWiny, cBlkx, cBlky, cCelx, cCely, cSrdx, cSrdy)
descriptors = zeros((round((cWinx-cBlkx)/cSrdx)+1)*(round((cWiny-cBlky)/cSrdy)+1)*(cBlkx/cCellx)*(cBlky/cCelly),9);
[m,n] = size(Grad);
%PixData: 信息数据
%(1,:)  histOfs[0]
%(2,:)  histOfs[1]
%(3,:)  histOfs[2]
%(4,:)  histOfs[3]
%(5,:)  histWeights[0]
%(6,:)  histWeights[1]
%(7,:)  histWeights[2]
%(8,:)  histWeights[3]
%(9,:)  gradOfs
%(10,:) qangleOfs
%(11,:) gradWeight
PixData = zeros(11,cBlkx*cBlky*3);
count1 = 0;
count2 = 0;
count4 = 0;
ncellx = floor(cBlkx/cCelx);
ncelly = floor(cBlky/cCely);
for i = 1:cBlkx
    for j = 1:cBlky
        cellX = (i+0.5)/cCelx-0.5;
        cellY = (j+0.5)/cCely-0.5;
        icellX0 = floor(cellX);
        icellY0 = floor(cellY);
        icellX1 = icellX0+1;
        icellY1 = icellY0+1;
        cellX = cellX-icellX0;
        cellY = cellY-icellY0;
        
        if ((icellX0 < ncellx) && (icellX1 < ncellx))
            if ((icellY0 < ncelly) && (icellY1 < ncelly))
                PixData(1,cBlkx*cBlky*2+count4) = (icellX0*ncelly+icellY0)*nbins;
                PixData(2,cBlkx*cBlky*2+count4) = (icellX1*ncelly+icellY0)*nbins;
                PixData(3,cBlkx*cBlky*2+count4) = (icellX0*ncelly+icellY1)*nbins;
                PixData(4,cBlkx*cBlky*2+count4) = (icellX1*ncelly+icellY1)*nbins;
                
                PixData(5,cBlkx*cBlky*2+count4) = (1.0-cellX)*(1.0-cellY);
                PixData(6,cBlkx*cBlky*2+count4) = (cellX)*(1.0-cellY);
                PixData(7,cBlkx*cBlky*2+count4) = (1.0-cellX)*(cellY);
                PixData(8,cBlkx*cBlky*2+count4) = (cellX)*(cellY);
                
                PixData(9,cBlkx*cBlky*2+count4) = (Grad(i,j)*i+j)*2;
                PixData(10,cBlkx*cBlky*2+count4) = (Angle(i,j)*i+j)*2;
                count4= count4+1;
            else
                if (icellY0 < ncelly)
                    icellY1 = icellY0;
                    cellY = 1.0 - cellY;
                end
                PixData(1,cBlkx*cBlky*1+count2) = (icellX0*ncelly+icellY1)*nbins;
                PixData(2,cBlkx*cBlky*1+count2) = (icellX1*ncelly+icellY1)*nbins;
                PixData(3,cBlkx*cBlky*1+count2) = 0;
                PixData(4,cBlkx*cBlky*1+count2) = 0;
                                      
                PixData(5,cBlkx*cBlky*1+count2) = (1.0-cellX)*(cellY);
                PixData(6,cBlkx*cBlky*1+count2) = (cellX)*(cellY);
                PixData(7,cBlkx*cBlky*1+count2) = 0;
                PixData(8,cBlkx*cBlky*1+count2) = 0;
                
                PixData(9,cBlkx*cBlky*1+count2) = (Grad(i,j)*i+j)*2;
                PixData(10,cBlkx*cBlky*1+count2) = (Angle(i,j)*i+j)*2;
                count2= count2+1;
            end
        else
            if (icellX0 < ncellx)
                icellX1 = icellX0;
                cellX = 1.0 - cellX;
            end
            if ((icellY0 < ncelly) && (icellY1 < ncelly))
            	PixData(1,cBlkx*cBlky*1+count2) = (icellX1*ncelly+icellY0)*nbins;
                PixData(2,cBlkx*cBlky*1+count2) = (icellX1*ncelly+icellY1)*nbins;
                PixData(3,cBlkx*cBlky*1+count2) = 0;
                PixData(4,cBlkx*cBlky*1+count2) = 0;
                                      
                PixData(5,cBlkx*cBlky*1+count2) = (cellX)*(1.0-cellY);
                PixData(6,cBlkx*cBlky*1+count2) = (cellX)*(cellY);
                PixData(7,cBlkx*cBlky*1+count2) = 0;
                PixData(8,cBlkx*cBlky*1+count2) = 0;
                
                PixData(9,cBlkx*cBlky*1+count2) = (Grad(i,j)*i+j)*2;
                PixData(10,cBlkx*cBlky*1+count2) = (Angle(i,j)*i+j)*2;
                count2= count2+1;
            else
            	if (icellY0 < ncelly)
                    icellY1 = icellY0;
                    cellY = 1.0 - cellY;
                end
            	PixData(1,count1) = (icellX1*ncelly+icellY1)*nbins;
            	PixData(2,count1) = 0;
            	PixData(3,count1) = 0;
            	PixData(4,count1) = 0;
            	
            	PixData(5,count1) = (cellX)*(cellY);
            	PixData(6,count1) = 0;
            	PixData(7,count1) = 0;
            	PixData(8,count1) = 0;
            	
                PixData(9,count1) = (Grad(i,j)*i+j)*2;
                PixData(10,count1) = (Angle(i,j)*i+j)*2;
            	count1= count1+1;
            end
        end
    end
end
BlkData = zeros(2, ((cWinx-cBlkx)/cSrdx+1)*((cWiny-cBlky)/cSrdy+1));
if (count1+count2+count3 == cBlkx*cBlky)
    for j = 1:count2
    	%pixData(j+count1) = pixData(j+cBlkx*cBlky);	
    end
    for j = 1:count4
    	%pixData(j+count1+count2) = pixData(j+cBlkx*cBlky*2);	
    end
    count2 = count2+count1;
    count4 = count4+count2;
    
    for i = 1:((cWinx-cBlkx)/cSrdx+1)
    	for j = 1:((cWiny-cBlky)/cSrdy+1)
    	    BlkData(1,i*((cWiny-cBlky)/cSrdy+1)+j) = (i*((cWiny-cBlky)/cSrdy+1)+j)*(cCelx*cCely*nbins);
    	    BlkData(2,i*((cWiny-cBlky)/cSrdy+1)+j) = 
    	end
    end
else
end


end
%% ham直方图统计

%function descriptor = ham_hist(InImage, xstart, ystart, xoffset, yoffset)
function descriptors_ = ham_hist3(Grad,Angle, Window, Block, Cell, Stride, nbins)
x = 1;y = 2;
%get Grad size
[m,n] = size(Grad);
n = n/2;
GradSize(x) = m;
GradSize(y) = n;
% get Window size
[m,n] = size(Window);
WinSize(x) = m;
WinSize(y) = n;
% get block size
[m,n] = size(Block);
BlkSize(x) = m;
BlkSize(y) = n;
% get cell size 
[m,n] = size(Cell);
CellSize(x) = m;
CellSize(y) = n;
% get stride size 
[m,n] = size(Stride);
StrideSize(x) = m;
StrideSize(y) = n;

% the numbers of blocks
nblock(x) = (WinSize(x)-BlkSize(x))/StrideSize(x)+1;
nblock(y) = (WinSize(y)-BlkSize(y))/StrideSize(y)+1;
nblocktotal = nblock(x)*nblock(y);
% the numbers of cells
ncell(x) = floor(BlkSize(x)/CellSize(x));
ncell(y) = floor(BlkSize(y)/CellSize(y));
nblockhistsize = ncell(x)*ncell(y)*nbins;
% the weight of block
BlkWeight = zeros(BlkSize(x), BlkSize(y));
sigma = 0.5;
if (sigma >= 0)
    sigma = sigma;
else
    sigma = round(BlkSize(x)+BlkSize(y))/8;
end
scale = 1.0/(sigma*sigma*2);

dx = zeros(1,BlkSize(x));
dy= zeros(1,BlkSize(y));
dh = BlkSize(y)*0.5;
dw = BlkSize(x)*0.5;
for i = 1:BlkSize(x)
    dx(1,i) = i - dw;
    dx(1,i) = dx(1,i)*dx(1,i);
end
for i = 1:BlkSize(y)
    dy(1,i) = i - dh;
    dy(1,i) = dy(1,i)*dy(1,i);
end

for i = 1:BlkSize(x)
    for j = 1:BlkSize(y)
        BlkWeight(i,j) = exp(-(dx(1,i)+dy(1,j))*scale);
    end
end
%     // Initialize 2 lookup tables, pixData & blockData.
%     // Here is why:
%     //
%     // The detection algorithm runs in 4 nested loops (at each pyramid layer):
%     //  loop over the windows within the input image
%     //    loop over the blocks within each window
%     //      loop over the cells within each block
%     //        loop over the pixels in each cell
%     //
%     // As each of the loops runs over a 2-dimensional array,
%     // we could get 8(!) nested loops in total, which is very-very slow.
%     //
%     // To speed the things up, we do the following:
%     //   1. loop over windows is unrolled in the HOGDescriptor::{compute|detect} methods;
%     //         inside we compute the current search window using getWindow() method.
%     //         Yes, it involves some overhead (function call + couple of divisions),
%     //         but it's tiny in fact.
%     //   2. loop over the blocks is also unrolled. Inside we use pre-computed blockData[j]
%     //         to set up gradient and histogram pointers.
%     //   3. loops over cells and pixels in each cell are merged
%     //       (since there is no overlap between cells, each pixel in the block is processed once)
%     //      and also unrolled. Inside we use PixData[k] to access the gradient values and
%     //      update the histogram
%     //

%% 初始化PixData LookUP Table
PixData=zeros(11, BlkSize(x)*BlkSize(y)*3);
%PixData(1:4,:)为HistOfs，PixData(5:8,:)为WeightOfs，PixData(9,:)表示GradOfs
%PixData(10,:)为AngleOfs，PixData(11,:)为GradWeight
HistOfs0 = 1;HistOfs1 = 2;HistOfs2 = 3;HistOfs3 = 4;
HistWeight0 = 5;HistWeight0 = 6;HistWeight0 = 7;HistWeight0 = 8;
GradOfs=9;AngleOfs=10;GradWeight=11;

count1 = 0;
count2 = 0;
count4 = 0;

for i = 1:BlkSize(x)
    for j = 1:BlkSize(y)
        cellX = (i+0.5)/CellSize(x)-0.5;
        cellY = (j+0.5)/CellSize(y)-0.5;
        icellX0 = floor(cellX);
        icellY0 = floor(cellY);
        icellX1 = icellX0+1;
        icellY1 = icellY0+1;
        cellX = cellX-icellX0;
        cellY = cellY-icellY0;
        
        if ((icellX0 < CellSize(x)) && (icellX1 < CellSize(x)))
            if ((icellY0 < CellSize(y)) && (icellY1 < CellSize(y)))
                PixData(1,BlkSize(x)*BlkSize(y)*2+count4) = (icellX0*CellSize(y)+icellY0)*nbins;
                PixData(2,BlkSize(x)*BlkSize(y)*2+count4) = (icellX1*CellSize(y)+icellY0)*nbins;
                PixData(3,BlkSize(x)*BlkSize(y)*2+count4) = (icellX0*CellSize(y)+icellY1)*nbins;
                PixData(4,BlkSize(x)*BlkSize(y)*2+count4) = (icellX1*CellSize(y)+icellY1)*nbins;
                
                PixData(5,BlkSize(x)*BlkSize(y)*2+count4) = (1.0-cellX)*(1.0-cellY);
                PixData(6,BlkSize(x)*BlkSize(y)*2+count4) = (cellX)*(1.0-cellY);
                PixData(7,BlkSize(x)*BlkSize(y)*2+count4) = (1.0-cellX)*(cellY);
                PixData(8,BlkSize(x)*BlkSize(y)*2+count4) = (cellX)*(cellY);
                
                PixData(9,BlkSize(x)*BlkSize(y)*2+count4) = (Grad(i,j)*i+j)*2;
                PixData(10,BlkSize(x)*BlkSize(y)*2+count4) = (Angle(i,j)*i+j)*2;
                PixData(11,BlkSize(x)*BlkSize(y)*2+count4) = BlkWeight(i,j);
                count4= count4+1;
            else %if ((icellY0 < CellSize(y)) && (icellY1 < CellSize(y)))
                if (icellY0 < CellSize(y))
                    icellY1 = icellY0;
                    cellY = 1.0 - cellY;
                end
                PixData(1,BlkSize(x)*BlkSize(y)*1+count2) = (icellX0*CellSize(y)+icellY1)*nbins;
                PixData(2,BlkSize(x)*BlkSize(y)*1+count2) = (icellX1*CellSize(y)+icellY1)*nbins;
                PixData(3,BlkSize(x)*BlkSize(y)*1+count2) = 0;
                PixData(4,BlkSize(x)*BlkSize(y)*1+count2) = 0;
                                      
                PixData(5,BlkSize(x)*BlkSize(y)*1+count2) = (1.0-cellX)*(cellY);
                PixData(6,BlkSize(x)*BlkSize(y)*1+count2) = (cellX)*(cellY);
                PixData(7,BlkSize(x)*BlkSize(y)*1+count2) = 0;
                PixData(8,BlkSize(x)*BlkSize(y)*1+count2) = 0;
                
                PixData(9,BlkSize(x)*BlkSize(y)*1+count2) = (Grad(i,j)*i+j)*2;
                PixData(10,BlkSize(x)*BlkSize(y)*1+count2) = (Angle(i,j)*i+j)*2;
                PixData(11,BlkSize(x)*BlkSize(y)*1+count2) = BlkWeight(i,j);
                count2= count2+1;
            end %if ((icellY0 < CellSize(y)) && (icellY1 < CellSize(y)))
        else %if ((icellX0 < CellSize(x)) && (icellX1 < CellSize(x)))
            if (icellX0 < CellSize(x))
                icellX1 = icellX0;
                cellX = 1.0 - cellX;
            end
            if ((icellY0 < ncelly) && (icellY1 < ncelly))
            	PixData(1,BlkSize(x)*BlkSize(y)*1+count2) = (icellX1*CellSize(y)+icellY0)*nbins;
                PixData(2,BlkSize(x)*BlkSize(y)*1+count2) = (icellX1*CellSize(y)+icellY1)*nbins;
                PixData(3,BlkSize(x)*BlkSize(y)*1+count2) = 0;
                PixData(4,BlkSize(x)*BlkSize(y)*1+count2) = 0;
                                      
                PixData(5,BlkSize(x)*BlkSize(y)*1+count2) = (cellX)*(1.0-cellY);
                PixData(6,BlkSize(x)*BlkSize(y)*1+count2) = (cellX)*(cellY);
                PixData(7,BlkSize(x)*BlkSize(y)*1+count2) = 0;
                PixData(8,BlkSize(x)*BlkSize(y)*1+count2) = 0;
                
                PixData(9,BlkSize(x)*BlkSize(y)*1+count2) = (Grad(i,j)*i+j)*2;
                PixData(10,BlkSize(x)*BlkSize(y)*1+count2) = (Angle(i,j)*i+j)*2;
                PixData(11,count2) = BlkWeight(i,j);
                count2= count2+1;
            else %if ((icellY0 < ncelly) && (icellY1 < ncelly))
            	if (icellY0 < CellSize(y))
                    icellY1 = icellY0;
                    cellY = 1.0 - cellY;
                end
            	PixData(1,count1) = (icellX1*CellSize(y)+icellY1)*nbins;
            	PixData(2,count1) = 0;
            	PixData(3,count1) = 0;
            	PixData(4,count1) = 0;
            	
            	PixData(5,count1) = (cellX)*(cellY);
            	PixData(6,count1) = 0;
            	PixData(7,count1) = 0;
            	PixData(8,count1) = 0;
            	
                PixData(9,count1) = (Grad(i,j)*i+j)*2;
                PixData(10,count1) = (Angle(i,j)*i+j)*2;                
                PixData(11,count1) = BlkWeight(i,j);
            	count1= count1+1;
            end %if ((icellY0 < ncelly) && (icellY1 < ncelly))
        end %if ((icellX0 < CellSize(x)) && (icellX1 < CellSize(x)))
    end %end of for j = 1:BlkSize(y)
end %end of i = 1:BlkSize(x)

if (count1+count2+count4 == BlkSize(x)*BlkSize(y))
    for j = 1:count2
        for i =1:11
            PixData(i,j+count1) = PixData(i,j+BlkSize(x)*BlkSize(y));	
        end
    end
    for j = 1:count4
        for i = 1:11
            PixData(i,j+count1+count2) = PixData(i,j+BlkSize(x)*BlkSize(y)*2);	
        end
    end
    count2 = count2+count1;
    count4 = count4+count2;

%% 初始化BlkData LookUP Table   
    histOfs = 1;
    imgOffsetx = 2;
    imgOffsety = 3;
    BlkData = zeros(3, nblocktotal);
    for i = 1:nblock(x)
        for j = 1:nblock(y)
            %块 的直方图偏移量
            BlkData(histOfs, (i-1)*nblock(y)+j) = ((i-1)*nblock(y)+j-1)*nblockhistsize;
            %块 的位置偏移量
            BlkData(imgOffsetx, (i-1)*nblock(y)+j) = (i-1)*StrideSize(x)+1;
            BlkData(imgOffsety, (i-1)*nblock(y)+j) = (j-1)*StrideSize(y)+1;
        end
    end
%% Calculate Histgram in Window
% getWindow block numbers
rho = 0;
s = rho;
for i = 1:nblock(x)
    for j = 1:nblock(y)
        Blockpoint(x) = BlkData(2,(i-1)*nblock(y)+j);
        Blockpoint(y) = BlkData(3,(i-1)*nblock(y)+j);
        hist = getBlock(Blockpoint, Blockdata);
        for k = 1:4:nblockhistsize
            s = s + hist(k)*svmVec(k)+hist(k+1)*svmVec(k+1)+hist(k+2)*svmVec(k+2)+hist(k+3)*svmVec(k+3);
        end
        for k = 1:nblockhistsize
            s = s + hist(k)*svmVec(k);
        end
    end %for j = 1:nblock(y)
end %for i = 1:nblock(x)
if (s >= hitThreshold)
    
end

else %if (count1+count2+count3 == BlkSize(x)*BlkSize(y))
    %Nothing todo
end %if (count1+count2+count3 == BlkSize(x)*BlkSize(y))
%%
% getBlock function
function blockHist=getBlock(BlockPoint, BlockData)
    BlockHist = zeros(4, nblockhistsize);
% for( k = 0; k < C1; k++ )
% {
%     const PixData& pk = _pixData[k];
%     const float* const a = gradPtr + pk.gradOfs;
%     float w = pk.gradWeight*pk.histWeights[0];
%     const uchar* h = qanglePtr + pk.qangleOfs;
%     int h0 = h[0], h1 = h[1];
% 
%     float* hist = blockHist + pk.histOfs[0];
%     float t0 = hist[h0] + a[0]*w;
%     float t1 = hist[h1] + a[1]*w;
%     hist[h0] = t0; hist[h1] = t1;
% }
%PixData(1:4,:)为HistOfs，PixData(5:8,:)为WeightOfs，PixData(9,:)表示GradOfs
%PixData(10,:)为AngleOfs，PixData(11,:)为GradWeight
    for k = 1:count1
        
    end %for k = 1:count1
    for k = 1:count2
        
    end %for k = 1:count2
    for k = 1:count4
        
    end %for k = 1:count4
end
end
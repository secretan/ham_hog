%% 使用Sobel算法进行梯度（）的大小和方向进行统计
% OutImage大小统计
% Vetor梯度方向
function [OutImage, Vetor] = my_sobel(InImage)
[m,n] = size(InImage);
parray = zeros(1,9);
paramx = [-1 0 1 -2 0 2 -1 0 1];
paramy = [1 2 1 0 0 0 -1 -2 -1];
%存放大小值
G = zeros(m,n);
%存放梯度方向
O = zeros(m,n);
%% 纵向处理与横向处理
for i = 1:m
    for j = 1:n
        if i == 1 
            if j == 1
                parray(1,1) = InImage(i,j);
                parray(1,2) = InImage(i,j);
                parray(1,3) = InImage(i,j);
                parray(2,1) = InImage(i,j);
                parray(3,1) = InImage(i,j);

                parray(3,2) = InImage(i,j+1);
                parray(3,3) = InImage(i+1,j+1);
            elseif j == n
                parray(1,1) = InImage(i,j);
                parray(2,1) = InImage(i,j);
                parray(3,1) = InImage(i,j);
                parray(3,3) = InImage(i,j);
                parray(3,2) = InImage(i,j);

                parray(1,2) = InImage(i,j-1);
                parray(1,3) = InImage(i+1,j-1);
            else                    
                parray(1,1) = InImage(i,j-1);
                parray(2,1) = InImage(i,j);
                parray(3,1) = InImage(i,j+1);

                parray(1,2) = InImage(i,j-1);
                parray(1,3) = InImage(i+1,j-1);
                parray(3,3) = InImage(i+1,j+1);
                parray(3,2) = InImage(i,j+1);
            end
            parray(2,2) = InImage(i,j);
            parray(2,3) = InImage(i+1,j);
        elseif i==m
            if j == 1
                parray(1,1) = InImage(i,j);
                parray(1,2) = InImage(i,j);
                parray(1,3) = InImage(i,j);
                parray(2,3) = InImage(i,j);
                parray(3,3) = InImage(i,j);

                parray(3,2) = InImage(i,j+1);
                parray(3,1) = InImage(i-1,j+1);
            elseif j == n
                parray(1,3) = InImage(i,j);
                parray(2,3) = InImage(i,j);
                parray(3,3) = InImage(i,j);
                parray(3,1) = InImage(i,j);
                parray(3,2) = InImage(i,j);

                parray(1,1) = InImage(i-1,j-1);
                parray(1,2) = InImage(i,j-1);
            else                    
                parray(1,3) = InImage(i,j-1);
                parray(2,3) = InImage(i,j);
                parray(3,3) = InImage(i,j+1);

                parray(1,1) = InImage(i-1,j-1);
                parray(1,2) = InImage(i,j-1);
                parray(3,1) = InImage(i-1,j+1);
                parray(3,2) = InImage(i,j+1);
            end
            parray(2,2) = InImage(i,j);
            parray(2,1) = InImage(i-1,j);
        else
            if j == 1
                parray(1,1) = InImage(i-1,j);
                parray(1,2) = InImage(i,j);
                parray(1,3) = InImage(i+1,j);
                parray(3,1) = InImage(i-1,j+1);
                parray(3,2) = InImage(i,j+1);
                parray(3,3) = InImage(i+1,j+1);
            elseif j == n
                parray(3,1) = InImage(i-1,j);
                parray(3,2) = InImage(i,j);
                parray(3,3) = InImage(i+1,j);
                parray(1,1) = InImage(i-1,j-1);
                parray(1,2) = InImage(i,j-1);
                parray(1,3) = InImage(i+1,j-1);
            else                    
            parray(1,1) = InImage(i-1,j-1);
            parray(1,2) = InImage(i,j-1);
            parray(1,3) = InImage(i+1,j-1);
            parray(3,3) = InImage(i+1,j+1);
            parray(3,2) = InImage(i,j+1);
            parray(3,1) = InImage(i-1,j+1);
            end
            parray(2,1) = InImage(i-1,j);
            parray(2,2) = InImage(i,j);
            parray(2,3) = InImage(i+1,j);
%             parray(1,1) = InImage(i-1,j-1);
%             parray(2,1) = InImage(i-1,j);
%             parray(3,1) = InImage(i-1,j+1);
%             parray(1,2) = InImage(i,j-1);
%             parray(1,3) = InImage(i+1,j-1);
%             parray(3,3) = InImage(i+1,j+1);
%             parray(3,2) = InImage(i,j+1);
%             parray(2,2) = InImage(i,j);
%             parray(2,3) = InImage(i+1,j);
        end      
       %计算大小      
       %方法一
       tmpx = 0;
       tmpy = 0;
       for k = 1:9
           tmpx = tmpx+ paramx(1,k)*parray(1,k);
           tmpy = tmpy+ paramy(1,k)*parray(1,k);
       end     
       G(i,j) = abs(tmpx)+abs(tmpy);
       %方法二
       %G(i,j) = abs((parray(1,1)+2*parray(1,2)+parray(1,3))-(parray(1,7)+2*parray(1,8)+parray(1,9)))+abs((parray(1,3)+2*parray(1,6)+parray(1,9))-(parray(1,1)+2*parray(1,4)+parray(1,7)));
       %计算方向
       O(i,j) = atan(abs(tmpx/tmpy));
    end
end
OutImage = G;
Vetor = O;
end
function se2_VizP_SingleSubj(SunjNum)


% baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze'; %laptop
baseDir = '/Users/nkordjazi/Documents/SeqEye/SeqEye2/analyze';
load([baseDir , '/CMB_34_1.mat'])
CMB = CMB_34_1;
 load([baseDir , '/se2_TranProb.mat']);
N2 = input('Use the last IPI for 2nd/3rd order transitions? (y/n)' , 's');
switch N2
    case 'y'
        LastIPI = 1;
    otherwise
        LastIPI = 0;
end
 clear PoC2 PoC3 PoC4
 
t2_Nums_allh = t2_Nums_allh(SunjNum);
t3_Nums_allh = t3_Nums_allh(SunjNum );
t4_Nums_allh = t4_Nums_allh(SunjNum);
CMB = CMB_34_1;
%%%%%%%%%%%%%%%************************************* 1st order
for p1 = 1:5
    for p2 = 1:5
        i = find(ismember(CMB.comb2 , [p1 p2] , 'rows'));
        PoC2(p1 , p2) = t2_Nums_allh.All(i);
        for h = [1:8 13]
            MT2_C{h}(p1, p2) = t2_Nums(h).MeanChunked_IPI(i);
            MT2_R{h}(p1, p2) = t2_Nums(h).MeanRand_IPI(i);
            MT2_A{h}(p1, p2) = t2_Nums(h).MeanAll_IPI(i);
        end
    end
    PoC2(p1, :) = PoC2(p1, :)/sum(PoC2(p1, :)); % sets the sum of every row to 1
end
figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC2);
ylabel('Press 1')
xlabel('Press 2')
axis square
colorbar
set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
    'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
title('Probability of Occurence for fist-order transitions in All sequences')
for j = 1:length(ChunkNum2)
    y = CMB.comb2(ChunkNum2(j,1) , 1);
    x = CMB.comb2(ChunkNum2(j,1) , 2);
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    imagesc(MT2_C{h} , [200 600]);
    ylabel('Press 1')
    xlabel('Press 2')
    axis square
    set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
        'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
    title(['Average Movement time (Chunked) H = ' , num2str(h)])
    for j = 1:length(ChunkNum2)
        y = CMB.comb2(ChunkNum2(j,1) , 1);
        x = CMB.comb2(ChunkNum2(j,1) , 2);
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    colorbar
    
    imcount = imcount+1;
end

figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC2);
ylabel('Press 1')
xlabel('Press 2')
axis square
colorbar
set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
    'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
title('Probability of Occurence for fist-order transitions in All sequences')
for j = 1:length(ChunkNum2)
    y = CMB.comb2(ChunkNum2(j,1) , 1);
    x = CMB.comb2(ChunkNum2(j,1) , 2);
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    imagesc(MT2_R{h} , [200 600]);
    ylabel('Press 1')
    xlabel('Press 2')
    axis square
    set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
        'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
    title(['Average Movement time (Random) H = ' , num2str(h)])
    for j = 1:length(ChunkNum2)
        y = CMB.comb2(ChunkNum2(j,1) , 1);
        x = CMB.comb2(ChunkNum2(j,1) , 2);
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    colorbar
    imcount = imcount+1;
end


figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC2);
ylabel('Press 1')
xlabel('Press 2')
axis square
colorbar
set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
    'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
title('Probability of Occurence for fist-order transitions in All sequences')
for j = 1:length(ChunkNum2)
    y = CMB.comb2(ChunkNum2(j,1) , 1);
    x = CMB.comb2(ChunkNum2(j,1) , 2);
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    imagesc(MT2_A{h} , [200 600]);
    ylabel('Press 1')
    xlabel('Press 2')
    axis square
    set(gca , 'XTick'  , [1:5] , 'YTick' , [1:5] , 'XTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'} , ...
        'Box' , 'off' , 'YTickLabels' , {'Thumb' , 'Index' , 'Middle' , 'Forth' , 'Pinkie'}  , 'LineWidth', 0.001)
    title(['Average Movement time (All) H = ' , num2str(h)])
    for j = 1:length(ChunkNum2)
        y = CMB.comb2(ChunkNum2(j,1) , 1);
        x = CMB.comb2(ChunkNum2(j,1) , 2);
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    colorbar
    imcount = imcount+1;
end

%%%%%%%%%%%%%%%************************************* 2nd order

for p1 = 1:5
    for p23 = 1:25
        i = find(ismember(CMB.comb3 , [p1 CMB.comb2(p23 , :)] , 'rows'));
        PoC3(p1 , p23) = t3_Nums_allh.All(i);
        for h = [1:8 13]
            MT3_C{h}(p1, p23) = t3_Nums(h).MeanChunked_IPI(i);
            MT3_R{h}(p1, p23) = t3_Nums(h).MeanRand_IPI(i);
            MT3_A{h}(p1, p23) = t3_Nums(h).MeanAll_IPI(i);
        end
    end
    PoC3(p1, :) = PoC3(p1, :)/sum(PoC3(p1, :)); % sets the sum of every row to 1
end

figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC3);
ylabel('Press 1')
xlabel('Press 2 3')
%         axis square
colorbar
set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
    'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
    'XTickLabelRotation' , 45)
title('Probability of Occurence for 2nd-order transitions in All sequences')
for j = 1:length(ChunkNum3)
    y = CMB.comb3(ChunkNum3(j,1) , 1);
    x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
for j = 1:length(ChunkNum3_4)
    y = CMB.comb3(ChunkNum3_4(j,1) , 1);
    x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    if ~LastIPI
        imagesc(MT3_C{h} , [500 , 1300]); % for sum of IPIS
    else
        imagesc(MT3_C{h},[200 600]);  % for the last IPIS
    end
    ylabel('Press 1')
    xlabel('Press 2 3')
    %             axis square
    set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
        'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
        'XTickLabelRotation' , 45)
    title(['Average Movement time (Chunked) H = ' , num2str(h)])
    
    for j = 1:length(ChunkNum3)
        y = CMB.comb3(ChunkNum3(j,1) , 1);
        x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    for j = 1:length(ChunkNum3_4)
        y = CMB.comb3(ChunkNum3_4(j,1) , 1);
        x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
    end
    
    colorbar
    
    imcount = imcount+1;
end


figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC3);
ylabel('Press 1')
xlabel('Press 2 3')
%         axis square
colorbar
set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
    'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
    'XTickLabelRotation' , 45)
title('Probability of Occurence for 2nd-order transitions in All sequences')
for j = 1:length(ChunkNum3)
    y = CMB.comb3(ChunkNum3(j,1) , 1);
    x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
for j = 1:length(ChunkNum3_4)
    y = CMB.comb3(ChunkNum3_4(j,1) , 1);
    x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    if ~LastIPI
        imagesc(MT3_R{h} , [500 , 1300]);% for sum of IPIS
    else
        imagesc(MT3_R{h},[200 600]);  % for the last IPIS
    end
    ylabel('Press 1')
    xlabel('Press 2 3')
    %             axis square
    set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
        'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
        'XTickLabelRotation' , 45)
    title(['Average Movement time (Random) H = ' , num2str(h)])
    
    for j = 1:length(ChunkNum3)
        y = CMB.comb3(ChunkNum3(j,1) , 1);
        x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    for j = 1:length(ChunkNum3_4)
        y = CMB.comb3(ChunkNum3_4(j,1) , 1);
        x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
    end
    
    colorbar
    
    imcount = imcount+1;
end


figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC3);
ylabel('Press 1')
xlabel('Press 2 3')
%         axis square
colorbar
set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
    'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
    'XTickLabelRotation' , 45)
title('Probability of Occurence for 2nd-order transitions in All sequences')
for j = 1:length(ChunkNum3)
    y = CMB.comb3(ChunkNum3(j,1) , 1);
    x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
for j = 1:length(ChunkNum3_4)
    y = CMB.comb3(ChunkNum3_4(j,1) , 1);
    x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    if ~LastIPI
        imagesc(MT3_A{h} , [500 , 1300]);% for sum of IPIS
    else
        imagesc(MT3_A{h},[200 600]);  % for the last IPIS
    end
    ylabel('Press 1')
    xlabel('Press 2 3')
    %             axis square
    set(gca , 'XTick'  , [1:25] , 'YTick' , [1:5] , 'YTickLabels' , {'1' , '2' , '3' , '4' , '5'} , ...
        'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
        'XTickLabelRotation' , 45)
    title(['Average Movement time (All) H = ' , num2str(h)])
    
    for j = 1:length(ChunkNum3)
        y = CMB.comb3(ChunkNum3(j,1) , 1);
        x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3(j,1) , 2:3) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    for j = 1:length(ChunkNum3_4)
        y = CMB.comb3(ChunkNum3_4(j,1) , 1);
        x = find(ismember(CMB.comb2 , CMB.comb3(ChunkNum3_4(j,1) , 2:3) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'green')
    end
    
    colorbar
    
    imcount = imcount+1;
end

%%%%%%%%%%%%%%%************************************* 3rd order

for p12 = 1:25
    for p34 = 1:25
        i = find(ismember(CMB.comb4 , [CMB.comb2(p12 , :) CMB.comb2(p34 , :)] , 'rows'));
        PoC4(p12 , p34) = t4_Nums_allh.All(i);
        for h = [1:8 13]
            MT4_C{h}(p12, p34) = t4_Nums(h).MeanChunked_IPI(i);
            MT4_R{h}(p12, p34) = t4_Nums(h).MeanRand_IPI(i);
            MT4_A{h}(p12, p34) = t4_Nums(h).MeanAll_IPI(i);
        end
    end
    PoC4(p12, :) = PoC4(p12, :)/sum(PoC4(p12, :)); % sets the sum of every row to 1
end

figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC4);
ylabel('Press 1 2')
xlabel('Press 3 4')
axis square
colorbar
set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
    'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
    'XTickLabelRotation' , 45)
title('Probability of Occurence for 3rd-order transitions in All sequences')
for j = 1:length(ChunkNum4)
    y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
    x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    if ~LastIPI
        imagesc(MT4_C{h},[800 2700]);  % for sum of IPIS
    else
        imagesc(MT4_C{h},[200 600]);  % for the last IPIS
    end
    ylabel('Press1')
    xlabel('Press2')
    axis square
    set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' ,cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
        'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
        'XTickLabelRotation' , 45)
    title(['Average Movement time (Chunked) H = ' , num2str(h)])
    
    for j = 1:length(ChunkNum4)
        y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
        x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    
    colorbar
    
    imcount = imcount+1;
end

figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC4);
ylabel('Press 1 2')
xlabel('Press 3 4')
axis square
colorbar
set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
    'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
    'XTickLabelRotation' , 45)
title('Probability of Occurence for 3rd-order transitions in All sequences')
for j = 1:length(ChunkNum4)
    y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
    x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    if ~LastIPI
        imagesc(MT4_R{h},[800 2700]);% for sum of IPIS
    else
        imagesc(MT4_R{h},[200 600]);  % for the last IPIS
    end
    
    ylabel('Press1')
    xlabel('Press2')
    axis square
    set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' ,cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
        'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
        'XTickLabelRotation' , 45)
    title(['Average Movement time (Random) H = ' , num2str(h)])
    
    for j = 1:length(ChunkNum4)
        y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
        x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    
    colorbar
    
    imcount = imcount+1;
end


figure('color' , 'white')
subplot(4,3,2)
imagesc(PoC4);
ylabel('Press 1 2')
xlabel('Press 3 4')
axis square
colorbar
set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
    'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
    'XTickLabelRotation' , 45)
title('Probability of Occurence for 3rd-order transitions in All sequences')
for j = 1:length(ChunkNum4)
    y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
    x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
    rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
end
imcount = 4;
for h = [1:8 , 13]
    subplot(4,3,imcount)
    if ~LastIPI
        imagesc(MT4_A{h} , [800 2700]);% for sum of IPIS
    else
        imagesc(MT4_A{h},[200 600]);  % for the last IPIS
    end
    ylabel('Press1')
    xlabel('Press2')
    axis square
    set(gca , 'XTick'  , [1:25] , 'YTick' , [1:25] , 'YTickLabels' ,cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false) , ...
        'Box' , 'off' , 'XTickLabels' , cellfun(@num2str , num2cell(CMB.comb2,2),'UniformOutput', false)  , 'LineWidth', 0.001,...
        'XTickLabelRotation' , 45)
    title(['Average Movement time (All) H = ' , num2str(h)])
    
    for j = 1:length(ChunkNum4)
        y = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 1:2) , 'rows'));
        x = find(ismember(CMB.comb2 , CMB.comb4(ChunkNum4(j,1) , 3:4) , 'rows'));
        rectangle('Position',[x - .5,y - .5,1,1],'LineWidth',2,'LineStyle',':' , 'EdgeColor' , 'red')
    end
    
    colorbar
    
    imcount = imcount+1;
end

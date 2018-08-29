%-------------------------------------------------------------------------%
% 
% This is a script to test the consistency of the noise model of the
% function data3Dgenerator by comparing the average value of each pixel
% with its variance.
% 
% !! ATTENTION1: if you want to look at a frame (of the simulated dataset 
%               stack) different from the first one you have to modify the
%               numFr parameter in the data3DGeneretor function and put it
%               equal to zero. 
% !! ATTENTION2: in the function data3Dgenerator comment the generator of
%                gaussian noise and then put                      
%                groundTruth(:,2)==disiredFrame instead of 0 
% !! ATTENTION3: remember to comment with % 
%                line 25:  folder_name = uigetdir ..
%                line 41:  figurePrint(art ..
%-------------------------------------------------------------------------%
nFr=100;
row=256;
col=256;
artificialFrame1 = zeros(row, col, nFr);
for i=1:100
    %groundTest = groundTruth;
    %groundTest(groundTest(:,2)~= 1, :)=[];
    artificialFrame1(:,:,i) = data3DGenerator( groundTruth(groundTruth(:, 2)==0, :), 256, 256, 80, 3.9, 500,  1);
end
meanFr1= mean(artificialFrame1,3);
stdFr1= std(artificialFrame1, 0,3);

figure, imagesc(stdFr1)
title(['Standard deviation of over ' num2str(nFr ) ' sample']);
figure, imagesc(meanFr1)
title(['Mean over ' num2str(nFr) ' sample']);
figure, imagesc(stdFr1.^2)
title(['Variance over ' num2str(nFr) ' sample']);
figure, imagesc(stdFr1.^2./meanFr1)
title(['Variance/Mean over ' num2str(nFr) ' sample']);
%-------------------------------------------------------------------------%
%-            THIS IS THE MAIN FUNCTION: 'artificialDHDataGenerator'     -%
%-------------------------------------------------------------------------%
% This function is a tool for generating a STORM Double Helix dataset 
% ( images of single molecules (SMs) ).
%
% This function gets in input:
%
% - calibrationDataPath: a stack of tiff files with the images of a bead at 
%   different z positions
% - zVect: a vector which encode the z position of each image in the 
%          calibration stack
%
% - groundTruthMatrix: a set of (x, y, z) coordinates
% - calibrationDarkFramePath
% - groundDarkFramePath
% - M: magnification factor
% - row: how many row the artifical images must have  
% - col: how many column the artificial images must have
%
% The function returns: 
%
% - a stack of tiff images

function [artificialDHdataMatrix] = artificialDHDataGenerator(calibrationDataPath, zVect, xyMatrix, groundTruthMatrix, calibrationDarkFramePath, groundDarkFramePath, M, row, col)

    % Loads image stacks of calibration images
    [rawTemplImagesStack, imgInfo, num_frames] = imageLoad(calibrationDataPath);
    [darkTemplImagesStack, imgInfo, num_frames] = imageLoad(calibrationDarkFramePath);
    [groundDarkStack, imgInfo, num_frames] = imageLoad(groundDarkFramePath);
    
    % Check inputs
    checkInput(row, col, max(groundTruthMatrix, 3), max(groundTruthMatrix, 4), M);
    
    % Templates dark count subtraction 
    imagesStack = darkSubtractor(rawTemplImagesStack, darkTemplImagesStack);
    
    % Generates the dataset and adds a source of noise
    artificialDHdataMatrix = dataGenerator(imagesStack, zVect, xyMatrix, groundTruthMatrix, M, row, col, groundDarkStack );
    
    % Saves on a tiff file the stack of artificial STORM DH data
    figurePrint(artificialDHdataMatrix);
end

%-------------------------------------------------------------------------%
%-                      THESE ARE THE OTHER METHODS                      -%
%-------------------------------------------------------------------------%

%__________________________________________________________________________
function figurePrint(imagesStack)

    for k = 1: size(imagesStack, 3)
    % Print the artificial data stack
        imwrite(uint16(round(imagesStack(:, :, k))), ['C:\Users\aarchett\Documents\MATLAB\3dDH\ArtificialAnnaWithNoiseImages' datestr(now,'yyyymmdd HHMM') '.tif'], 'writemode', 'append');
    end
end
 

%__________________________________________________________________________
function [images, imgInfo, num_frames] = imageLoad(fName)
    
    % For each path loads an images stack file 'fName' and return a 3D 
    % matrix contening the images pixel values in the firts two dimentions
    % (x and y) and a different image in each z value.
    
    % Grab file info.
    imgInfo = imfinfo(fName);
    num_frames = numel(imgInfo);
    
    % Give image width and heigh
    colNum = imgInfo(1).Width;
    rowNum = imgInfo(1).Height;
    
    images = zeros(rowNum, colNum,  num_frames, 'uint16');
    
    for k = 1:num_frames
   
    % Actually load images.    
        images(:,:,k) = imread(fName, k);
    end
end

%__________________________________________________________________________
function imagesStack = darkSubtractor(rawImagesStack, darckImagesStack)
    
    averagedDark = uint16(mean(darckImagesStack, 3));
    imagesStack = rawImagesStack - repmat(averagedDark, [1, 1, size(rawImagesStack, 3)]) ;

end

%__________________________________________________________________________
function [dataset] = dataGenerator(templateStack, zVect, xyMatrix, groundTruthMatrix, M, rows, cols, groundDarkStack)
    
    % Preallocate matrix
    dataset = uint16( zeros( rows, cols, max(groundTruthMatrix(:, 2)) +1 ) );
    %dataset = uint16(zeros( round(max(groundTruthMatrix(:, 4)*M) + 10), round(max(groundTruthMatrix(:, 3)*M) + 10), max(groundTruthMatrix(:, 2)) ));
    
    % For each molecule in the groundTruthMatrix
    for molecIdx = 1: size(groundTruthMatrix, 1)

        % Grab the x, y, z position of the molecule per a magnification factor M
        zPos = round(groundTruthMatrix(molecIdx, 5)*M);
        xPos = round(groundTruthMatrix(molecIdx, 3)*M);
        yPos = round(groundTruthMatrix(molecIdx, 4)*M);
        framePos = groundTruthMatrix(molecIdx, 2) + 1;
        
        zMax = max(groundTruthMatrix(:, 5)*M);
        
        % Ridefine the z position with the z calibration vector
        newZValue = round( zPos*zVect(length(zVect))/zMax );
        
        % Save the position of the corresponding template
        % This finds the value in zVect which is closest to the newZPos value I am calling.
        [value, templIdx] = min(abs(zVect-newZValue));
        
        %templIdx = find(zVect == newZPos);
        XtemplCenter = xyMatrix(templIdx, 1) ;
        YtemplCenter = xyMatrix(templIdx, 2);
        templHeight = size( templateStack(:, :, templIdx), 1);
        templWidth = size( templateStack(:, :, templIdx), 2);
        
        % Assigne the template values to the dataset matrix at the molecule
        % position (making the template centre coincident with the molecule
        % position)
        dataset( yPos - YtemplCenter : templHeight + (yPos - YtemplCenter) -1, xPos - XtemplCenter : templWidth + (xPos - XtemplCenter)-1,  framePos) = templateStack(:, :, templIdx);
        
    end
    dataset = dataset + groundDarkStack./10 ;
end

%__________________________________________________________________________
function checkInput(row, col, xMax, yMax, M)
    
    if row < yMax*M
        error('The number of rows chosen is smaller than the maximum molecule y position');
    end
    if col < xMax*M
        error('The number of columns chosen is smaller than the maximum molecule x position');
    end
end


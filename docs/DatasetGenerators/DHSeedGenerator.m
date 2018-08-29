%-------------------------------------------------------------------------%
%-            THIS IS THE MAIN FUNCTION: 'DHSeedGenerator'               -%
%-------------------------------------------------------------------------%
% This function is a tool for generating a stack of seed. One seed for each
% z position.
%
% This function gets in input:
%
% - row: how many row the artifical images must have  
% - col: how many column the artificial images must have
% - zRange: the range in z [nm]
% - pixSize: the size of the camera pixel [nm]
% - I0: the number of photons emitted at the center of the molecule [photons]
% - sigma: the PSF [nm]
% - dist: distance between the splitted molecule images
% - ground: average background
% - numFr: the number of frame 
%
% The function returns: 
%
% - a stack of tiff images
%
% ATTENTION : the z axis is calibrated in this way:  zRange : 2Pi = z : angle
%             and starts with z=0 <----> angle=0

function [DHseedMatrix] = DHSeedGenerator( row, col, zRange, pixSize, I0, sigma, dist, ground, numFr)

    % Generates a base of background normally distributed
    base = poissrnd(ground, [row, col, numFr]);
    
    % Generates a groundTruth with 1 molecule per frame. All in the central
    % position
    groundTruthMatrix = zeros(numFr, 4);
    groundTruthMatrix(:, 2) = repmat( round(col/2)*pixSize, numFr, 1);
    groundTruthMatrix(:, 3) = repmat( round(row/2)*pixSize, numFr, 1);
    groundTruthMatrix(:, 4) = (1: zRange/numFr : zRange)';
    groundTruthMatrix(:, 1) = (1: 1: numFr)';
    
    % Trasforms the x,y molecule coordinate in the x,y DH coordinate
    [X1, X2, Y1, Y2, fr]  = coordTransformator(groundTruthMatrix, dist, zRange);
    
    % Generates the dataset and adds a source of noise
    DHseedMatrix = dataGenerator(X1, X2, Y1, Y2, fr, base, pixSize, I0, sigma, row, col);
    
    % Saves on a tiff file the stack of artificial STORM DH data
    figurePrint(DHseedMatrix);
end

%-------------------------------------------------------------------------%
%-                      THESE ARE THE OTHER METHODS                      -%
%-------------------------------------------------------------------------%

%__________________________________________________________________________
function figurePrint(imagesStack)

    for k = 1: size(imagesStack, 3)
    % Print the artificial data stack
        imwrite(uint16(round(imagesStack(:, :, k))), ['C:\Users\aarchett\Documents\MATLAB\DHseed' datestr(now,'yyyymmdd HHMM') '.tif'], 'writemode', 'append');
    end
end

%__________________________________________________________________________
function [X1, X2, Y1, Y2, fr] = coordTransformator(groundTruthMatrix, dist, zRange)

    % Computes the position of the two DH molecule images for each molecule
    % (in nm!!)
    X1 = groundTruthMatrix(:, 2) - (dist/2)*cos( 2*pi*groundTruthMatrix(:, 4)/zRange ) ;
    X2 = groundTruthMatrix(:, 2) + (dist/2)*cos( 2*pi*groundTruthMatrix(:, 4)/zRange ) ;
    Y1 = groundTruthMatrix(:, 3) - (dist/2)*sin( 2*pi*groundTruthMatrix(:, 4)/zRange ) ;
    Y2 = groundTruthMatrix(:, 3) + (dist/2)*sin( 2*pi*groundTruthMatrix(:, 4)/zRange ) ;
    fr = groundTruthMatrix(:, 1);
    
end

%__________________________________________________________________________
function [base] = dataGenerator(X1, X2, Y1, Y2, fr, base, pixSize, I0, sigma, row, col)
    
    % For each molecule in the groundTruthMatrix 
    for molIdx = 1 : length(X1) 
        
        % Finds the coordinates of the pixel closest to the molecule position
        % [in pixel]
        X1p = round(X1(molIdx) / pixSize);
        Y1p = round(Y1(molIdx) / pixSize);
        
        X2p = round(X2(molIdx) / pixSize);
        Y2p = round(Y2(molIdx) / pixSize);
        
        frm = fr(molIdx);
        
        for rowId = Y1p - row/4 : Y1p + row/4
            for colId = X1p - col/4 : X1p + col/4
                
            % Computes the coordinates of the pixel position in the molecule
            % center frame reference [in nm]
            X1d = abs(X1(molIdx) - colId*pixSize);
            Y1d = abs(Y1(molIdx) - rowId*pixSize);

            % In each pixel of the simulated image there is a number of photons 
            % Poisson distributed with an average given by the value of the Gaussian
            % PSF in that point
            Gauss1 = I0*exp(-(X1d^2 + Y1d^2)/(2*sigma^2)); 
            
            base(rowId, colId, frm) = base(rowId, colId, frm) + poissrnd(Gauss1);
            end 
        end   
        for rowId = Y2p - row/4 : Y2p + row/4
            for colId = X2p - col/4 : X2p + col/4
                
            % Computes the coordinates of the pixel position in the molecule
            % center frame reference [in nm]   
            X2d = abs(X2(molIdx) - colId*pixSize);
            Y2d = abs(Y2(molIdx) - rowId*pixSize);
    
            % In each pixel of the simulated image there is a number of photons 
            % Poisson distributed with an average given by the value of the Gaussian
            % PSF in that point
            Gauss2 = I0*exp(-(X2d^2 + Y2d^2)/(2*sigma^2));
            
            base(rowId, colId, frm) = base(rowId, colId, frm) + poissrnd(Gauss2);
            end 
        end
        
    end
end



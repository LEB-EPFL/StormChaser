%-------------------------------------------------------------------------%
%-            THIS IS THE MAIN FUNCTION: 'data3DGenerator'               -%
%-------------------------------------------------------------------------%
% This function is a tool for generating a STORM Double Helix dataset 
% ( images of single molecules (SMs) ).
%
% This function gets in input:
%
% - Mxy: magnification factor (multiplies the x, y position in pix by a factor of Mxy)
% - Mz: magnification factor (multiplies the z position in pix by a factor of Mz)
% - zShift: shift the z coordinates by a constant value zShift only for
%   generating higher intensities (line 106) not all z!!!
% - groundTruthMatrix: a set of (n, fr, x[pix], y[pix], z[axial resolution]) coordinates (has to be a
%                      matlab matrix) fr = frame starting from 0
% - row: how many row the artifical images must have  
% - col: how many column the artificial images must have
% - pixSize: the size of the sample pixel [nm]
% - axialRes: axial resolution [nm] (4)
% - ground: average background 
% - PSF: enter 1 for a DHPSF, 2 for an ..?.. PSF
%
% The function returns: 
%
% - a stack of tiff images
%
% Input with the groundTruth given is: [artificialDHdataMatrix12] = data3DGenerator(2, 2, -1000, groundTruth, 512, 512, 120, 3.9, 100, 1);


function [artificialDataMatrix] = data3DGenerator(Mxy, Mz, zShift, groundTruthMatrix, row, col, pixSize, axialRes, ground, PSF)

    %load('Z:\Users\Anna-Archetti\DHproj\modelCheckImg\groundTruth.mat')
    %folder_name = uigetdir(matlabroot, 'Choose where to save the images with the artificial dataset');
    folder_name = 0;
    
    % Original Pixel size in the calibration bead [nm]!!
    origPixSize = 126; 
    
    % How big you want your PSF signal (number of sigma)
    sigmaNr = 6;
    
    % Grab the number of frame
    numFr = max(groundTruthMatrix(:, 2) + 1 );
    %  numFr = max(groundTruthMatrix(:, 2));
    
    % Generates a base of background normally distributed with mean =
    % ground and sigma = 10
      %base = normrnd(ground, 10, [row, col, numFr]);
      %base = zeros(row, col, numFr);
    
      
    % Generates a base of background poisson distributed with mean = ground
    base = poissrnd(ground, [row, col, numFr]);
    
    % Check if the number of row and col match with the max number of pixel 
    % in the groundTruthMatrix
    check(groundTruthMatrix, row, col, axialRes, Mz, pixSize, origPixSize, sigmaNr);
        
    % Generates the PSF and adds a source of poisson noise
    artificialDataMatrix = dataGenerator(Mxy, Mz, zShift, groundTruthMatrix, base, pixSize, axialRes, PSF, origPixSize, sigmaNr);
    
    % Saves on a tiff file the stack of artificial STORM DH data
    figurePrint(artificialDataMatrix, folder_name);
end

%-------------------------------------------------------------------------%
%-                      THESE ARE THE OTHER METHODS                      -%
%-------------------------------------------------------------------------%

%__________________________________________________________________________
function figurePrint(imagesStack, folder_name)

    for k = 1: size(imagesStack, 3)
    % Print the artificial data stack
        %imwrite(uint16(round(imagesStack(:, :, k))), [strcat( folder_name, '\artificialDatasetDH3Dbead', datestr(now,'yyyymmdd HHMM')) '.tif'], 'writemode', 'append');
        imwrite(uint16(round(imagesStack(:, :, k))), 'W:\LEB\Users\Anna-Archetti\DH-Challenge 2016\possible_psf-DH\DH_PSF_simulated\artificialPSF_16bit.tif', 'writemode', 'append');
        
        
        t = Tiff(['W:\LEB\Users\Anna-Archetti\DH-Challenge 2016\possible_psf-DH\DH_PSF_simulated\32bit_series\artificialPSF_32bit' num2str(k) '.tif'], 'w');
        tagstruct.ImageLength = size(imagesStack(:, :, k), 1);
        tagstruct.ImageWidth = size(imagesStack(:, :, k), 2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 32;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.RowsPerStrip    = 16;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.Software        = 'MATLAB';
        t.setTag(tagstruct);
        t.write( uint32( round(imagesStack(:, :, k)) ) );
        t.close();
        
        disp('Frame number:')
        disp(k)
        disp(' written')
    end
    
end

%__________________________________________________________________________
function [X1, X2, Y1, Y2, fr, z] = coordTransformator(Mxy, Mz, zShift, groundTruthMatrix, pixSize, axialRes, origPixSize)
    
    % Grab the number of molecules
    numMol = size(groundTruthMatrix,1);
    %z = groundTruthMatrix(:, 5).*axialRes.*Mz + repmat(zShift, numMol, 1);
    z = groundTruthMatrix(:, 5).*axialRes.*Mz ;
    
    % Define dist: distance between the splitted molecule images [pixel]
    %dist = 3;
    dist = (7E-07.*z.^2 - 0.0008.*z + 5.8).*(origPixSize/pixSize);
    
    % Computes the position of the two DH molecule images for each molecule
    % (in nm!!)
    % Angle in degree!!
    
    angle = -0.07.*z;
    % angle = -0.0687.*z + 5.2449;
    %angle = -0.085.*z + 171;
    
    
    % The pixels of a image are count from the top to the botton so angle
    % are swap respect to x axis when they will be "read"
    X1 = Mxy.*groundTruthMatrix(:, 3) - (dist./2).*cosd( angle ) ;
    X2 = Mxy.*groundTruthMatrix(:, 3) + (dist./2).*cosd( angle ) ;
    Y1 = Mxy.*groundTruthMatrix(:, 4) - (dist./2).*sind( angle ) ;
    Y2 = Mxy.*groundTruthMatrix(:, 4) + (dist./2).*sind( angle ) ;

    
    fr = groundTruthMatrix(:, 2) + 1;
   % fr = groundTruthMatrix(:, 2);
    
end

%__________________________________________________________________________
function [base] = DHPSFmodel(Mxy, Mz, zShift, groundTruthMatrix, base, pixSize, axialRes, origPixSize, sigmaNr)

    % Trasforms the x,y molecule coordinate in the x,y DH coordinates (everithing in rescaled pix)
    [X1, X2, Y1, Y2, fr, z]  = coordTransformator(Mxy, Mz, zShift, groundTruthMatrix, pixSize, axialRes, origPixSize);
    
    % Define sigma: the PSF [nm]
    %sigma = 1.25*pixSize;
    sigma1 = (4E-07.*z.^2 - 3E-05.*z + 1.45).*origPixSize;
    sigma2 = (9E-07*z.^2 - 0.0004.*z + 1.5).*origPixSize;
    
    % Define I0: the number of photons emitted at the center of the molecule [photons]
    I0 = -3.6.*(z + zShift) + 5000;
    
    %I0 = 10000./(2*pi*(sigma./pixSize).^2); 
   % I0 = repmat(400, 1, length(z));
   % I0 =-0.2*z+500;
   
    % For each molecule in the groundTruthMatrix     
    figure,
    for molIdx = 1 : length(X1) 
        
        % Finds the coordinates of the pixel closest to the molecule position
        % [in pixel] respect to the new image with pixel with a differepixSizent
        % size
        X1p = round(X1(molIdx));
        Y1p = round(Y1(molIdx));
        
        X2p = round(X2(molIdx));
        Y2p = round(Y2(molIdx));
        
        frm = fr(molIdx);
        
        PSFSigma = round( max(abs(4E-07.*z.^2 - 3E-05.*z + 1.45)*(origPixSize/pixSize)))*sigmaNr; % pixel 
        for rowId = Y1p - PSFSigma : Y1p + PSFSigma
            for colId = X1p - PSFSigma : X1p + PSFSigma
                
            % Computes the coordinates of the pixel position in the molecule
            % center frame reference [in nm]
            X1d = abs(X1(molIdx)*pixSize - colId*pixSize);
            Y1d = abs(Y1(molIdx)*pixSize - rowId*pixSize);

            % In each pixel of the simulated image there is a number of photons 
            % Poisson distributed with an average given by the value of the Gaussian
            % PSF in that point
            Gauss1 = ( I0(molIdx)/2 )*exp(-(X1d^2 + Y1d^2)/(2*sigma1(molIdx)^2)); 
            
            base(rowId, colId, frm) = base(rowId, colId, frm) + poissrnd(Gauss1);
            end 
        end   
        for rowId = Y2p - PSFSigma : Y2p + PSFSigma
            for colId = X2p - PSFSigma : X2p + PSFSigma
                
            % Computes the coordinates of the pixel position in the molecule
            % center frame reference [in nm]   
            X2d = abs(X2(molIdx)*pixSize - colId*pixSize);
            Y2d = abs(Y2(molIdx)*pixSize - rowId*pixSize);
    
            % In each pixel of the simulated image there is a number of photons 
            % Poisson distributed with an average given by the value of the Gaussian
            % PSF in that point
            Gauss2 = ( I0(molIdx)/2 )*exp(-(X2d^2 + Y2d^2)/(2*sigma2(molIdx)^2));
            
            base(rowId, colId, frm) = base(rowId, colId, frm) + poissrnd(Gauss2);
            end
        end
        % DEBUG
        imagesc(base(:, :, frm));
        drawnow;
    end

end

%__________________________________________________________________________
function [dataset] = dataGenerator(Mxy, Mz, zShift, groundTruthMatrix, base, pixSize, axialRes, PSF, origPixSize, sigmaNr)  
    if PSF == 1 
        [dataset]= DHPSFmodel(Mxy, Mz, zShift, groundTruthMatrix, base, pixSize, axialRes, origPixSize, sigmaNr);
    else
        error('You have to implement another PSF :P!')
    end
end

%__________________________________________________________________________
function check(groundTruthMatrix, row, col, axialRes, Mz, pixSize, origPixSize, sigmaNr)
    
    xMax = max(groundTruthMatrix(:, 3));
    yMax = max(groundTruthMatrix(:, 4));
    
    z = groundTruthMatrix(:, 5).*axialRes.*Mz ;
    sigma1Max = max(abs(4E-07.*z.^2 - 3E-05.*z + 1.45))*(origPixSize/pixSize);
    sigma2Max = max(abs(9E-07*z.^2 - 0.0004.*z + 1.5))*(origPixSize/pixSize);
    
    distMax = max(abs(7E-07.*z.^2 - 0.0008.*z + 5.8))*(origPixSize/pixSize);
    if  col < xMax 
        error('The numer of columns chosen is smaller than the maximun x value');
    end
    if col < sigma1Max*sigmaNr && col < sigma2Max*sigmaNr
        error(['The number of columns chosen is smaller than the double sigma of the PSF that is:' num2str(round( max(max(sigma1), max(sigma2))*2 )) ])
    end 
    if  row < yMax
        error('The numer of rows chosen is smaller than the maximun y value');
    end
    if row < sigma1Max*sigmaNr && row < sigma2Max*sigmaNr
        error(['The number of rows chosen is smaller than the sigma of the PSF that is:' num2str(round( max(max(sigma1), max(sigma2))*2 ))])
    end 
    if row < distMax*2
        error(['The number of rows chosen is smaller than the distance of the two lobes of the PSF that is:' num2str(round( max(dist))*2 )])
    end
    if col < distMax*2
        error(['The number of columns chosen is smaller than the distance of the two lobes of the PSF that is:' num2str(round( max(dist))*2 )])
    end
    if row < distMax*2 + max(sigma1Max*sigmaNr, sigma2Max*sigmaNr)
         error(['The number of rows chosen is smaller than the distance of the two lobes of the PSF + the sigma of the PSF, that is, smaller than:' num2str(round( distMax*2 + max(sigma1Max*2, sigma2Max*2) ))])
    end
    if col < distMax*2 + max(sigma1Max*sigmaNr, sigma2Max*sigmaNr)
         error(['The number of columns chosen is smaller than the distance of the two lobes of the PSF + the sigma of the PSF, that is, smaller than:' num2str(round( distMax*2 + max(sigma1Max*2, sigma2Max*2) ))])
    end
    
    
end



function [y1, y2, y3] = enf(x, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency)

%enf.m Summary of this function goes here

%First segment the input 'x' into blocks using BlockSize and Overlap
%script will use audioread to create vector x and sends through function
increment = BlockSize - (BlockSize*Overlap);

row = [1:BlockSize];
col = [0:increment:((length(x))-(BlockSize+1))];
[X,Y] = meshgrid(row,col);

indexMatrix = X + Y;
inputSegment = x(indexMatrix); %segmented matrix
[rowLen,~] = size(inputSegment);

%Need to window, zeropad, and perform fft on each segment:
%Create a windowing function and apply to each frame note window returns
%column vector and segmented input is by rows
windowRow = transpose(Window);

%Repmat(matrix, row multiplier, col multiplier) - want window to be same
%size as segmented matrix
d = repmat(windowRow,rowLen,1);

%Element wise multiplication between the segmented input and repmatted window
windowedFrame = inputSegment.*d;

%Zero pad the frame by creating array of zeros and concatenate. Matlab cat
%function only works with arrays of same size
zeroos = zeros(1, ZeroPad);
zeroArray = repmat(zeroos,rowLen,1);
afterCat = cat(2, windowedFrame, zeroArray);

%Fianlly find fft on each block
final = fft(afterCat,length(afterCat),2);

%Now need to get frequencies of interest or else compooter go poof
%Consider new block size after zero padding

N = BlockSize + ZeroPad;
frequencies = Fs*[0:N-1]/N;

%k is the indices of the desired frequencies
k1 = Frequency-1;
k2 = Frequency+1;
k = find(frequencies>=(k1) & frequencies <=(k2));
newVar = frequencies(k);

%Getting y1, the surface plot (magnitude)
surfPlot= abs(final(:,k)); 
y1 = surfPlot;

%Getting y2, the maximum magnitude
[~,I] = max(surfPlot,[],2);
y2 = newVar(I);

%Getting y3, the weighted magnitude 
[row, col] = size(surfPlot);
desiredFreq = repmat(k,row,1);
numerator = sum((frequencies(desiredFreq).*surfPlot),2);
denominator = sum(surfPlot,2);
weightedMag = numerator./denominator;
y3 = weightedMag;

end


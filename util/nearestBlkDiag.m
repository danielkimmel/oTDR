
function [nBlkDiagMtx, blkIxs, blkSizes, blkValues] = nearestBlkDiag(Mtx)
N = length(Mtx);
Er = norm(Mtx - diag(diag(Mtx)), 'fro');
if N>1
    for blkSize = 2:N
        numPatterns = N- blkSize+1;
        
        for i = 1:numPatterns
            blkMtx = diag(diag(Mtx));
            blkMtx_i = Mtx(i:i+blkSize-1, i:i+blkSize-1);
            blkVal = mean(blkMtx_i(:));
            blkMtx(i:i+blkSize-1, i:i+blkSize-1) = blkVal;
            Er_blk(i) = norm(Mtx - blkMtx, 'fro');      
        end
        
        Er(blkSize) = min(Er_blk);
    end
    [~, blkSize] = min(Er);
    numPatterns = N- blkSize+1; 
    for i = 1:numPatterns
        blkMtx = diag(diag(Mtx));
        blkMtx_i = Mtx(i:i+blkSize-1, i:i+blkSize-1);
        blkVal = mean(blkMtx_i(:));
        blkMtx(i:i+blkSize-1, i:i+blkSize-1) = blkVal;
        Er_blk(i) = norm(Mtx - blkMtx, 'fro');      
    end
    [~, blkIx] = min(Er_blk);
else
     blkIx = 1;
     blkSize = 1;
end
[resMtx1, resMtx2] = getResidualMtx(Mtx, blkIx, blkSize);
if ~isempty(resMtx1)
    [~, blkIx1, blkSize1] = nearestBlkDiag(resMtx1);
    blkIx1 = [blkIx1(:); blkIx];
    blkSize1 = [blkSize1(:); blkSize];
else
    blkIx1 = blkIx;
    blkSize1 = blkSize;
end

if ~isempty(resMtx2)
    [~, blkIx2, blkSize2] = nearestBlkDiag(resMtx2);
    blkIxs = [blkIx1(:); blkIx2(:)+blkIx+blkSize-1];
    blkSizes = [blkSize1(:); blkSize2(:)];
else
    blkIxs = blkIx1(:);
    blkSizes = blkSize1(:);
end

blkIxs = blkIxs(:);
blkSizes = blkSizes(:);
nBlkDiagMtx = zeros(N);
blkValues = nan(length(blkIxs),1);
for i = 1:length(blkIxs)
    Mtx_i = Mtx(blkIxs(i):blkIxs(i)+blkSizes(i)-1, blkIxs(i):blkIxs(i)+blkSizes(i)-1);
    blkValues(i) = mean(Mtx_i(:));
    nBlkDiagMtx(blkIxs(i):blkIxs(i)+blkSizes(i)-1, blkIxs(i):blkIxs(i)+blkSizes(i)-1) = blkValues(i); 
    blkValues(i) = mean(mean(Mtx_i-diag(diag(Mtx_i))));
end
end

    

function [resMtx1, resMtx2] = getResidualMtx(Mtx, blkIx, blkSize)
    N = length(Mtx);
    resMtx1 = Mtx(1:blkIx-1, 1:blkIx-1);
    
    resMtx2 = Mtx(blkIx+blkSize:N, blkIx+blkSize:N);
end


            
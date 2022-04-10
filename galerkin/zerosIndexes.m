function retVal = zerosIndexes(vector)
    zci = @(v) find(diff(sign(v)));

    neighborhood = 3;

    diffs = diff(vector);
    vectorZeroInds = zci(vector);

    for i = 1:1:size(vectorZeroInds, 2)
        diffMinInd = vectorZeroInds(i) - neighborhood;

        if (diffMinInd < 1)
            diffMinInd = 1;
        end

        diffMaxInd = vectorZeroInds(i) + neighborhood;

        if (diffMaxInd > size(diffs, 2))
            diffMaxInd = size(diffs, 2);
        end

        diffAbsInterval = abs(diffs(diffMinInd:diffMaxInd));
        diffAbsIntervalMean = mean((diffAbsInterval), 2, 'omitnan');

        numberOfEjections = size(diffAbsInterval(diffAbsInterval < diffAbsIntervalMean / 1000), 2) + ...
            size(diffAbsInterval(diffAbsInterval > diffAbsIntervalMean * 1000), 2);
        numberOfNans = size(diffAbsInterval(isnan(diffAbsInterval)), 2);

        if (numberOfNans > 0 || numberOfEjections > 0)
            vectorZeroInds(i) = NaN;
        end

    end

    retVal = vectorZeroInds(~isnan(vectorZeroInds));

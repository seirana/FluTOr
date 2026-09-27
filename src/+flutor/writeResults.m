function paths = writeResults(outputDirectory, runName, output, metadata)
%WRITERESULTS Persist FluTOr results after computation completes.

    outputDirectory = string(outputDirectory);
    runName = string(runName);
    if strlength(runName) == 0
        error("flutor:InvalidRunName", "runName must not be empty.");
    end

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    stem = regexprep(runName, "[^A-Za-z0-9._-]+", "_");
    tradeoffPath = fullfile(outputDirectory, stem + "_tradeoffs.csv");
    resultPath = fullfile(outputDirectory, stem + "_analysis.mat");
    metadataPath = fullfile(outputDirectory, stem + "_metadata.json");

    writetable(output.tradeoffResult.tradeoffs, tradeoffPath);

    analysis = output; %#ok<NASGU>
    save(resultPath, "analysis", "-v7.3");

    metadata.solverRuns = output.tradeoffResult.solverRuns;
    metadata.nTradeoffs = height(output.tradeoffResult.tradeoffs);
    metadata.preprocessing = output.preprocessing;
    metadata.couplingDiagnostics = output.coupling.diagnostics;

    fid = fopen(metadataPath, "w");
    if fid < 0
        error("flutor:OutputWriteFailed", "Could not open %s for writing.", metadataPath);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, "%s", jsonencode(metadata, "PrettyPrint", true));
    clear cleanup;

    paths = struct( ...
        "tradeoffs", string(tradeoffPath), ...
        "analysis", string(resultPath), ...
        "metadata", string(metadataPath));
end

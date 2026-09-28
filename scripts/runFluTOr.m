function output = runFluTOr(modelFile, biomassReaction, outputDirectory, boundsFile)
%RUNFLUTOR File-driven entry point for the maintained FluTOr pipeline.

    arguments
        modelFile (1, 1) string
        biomassReaction (1, 1) string
        outputDirectory (1, 1) string = "artifacts"
        boundsFile (1, 1) string = ""
    end

    repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
    addpath(fullfile(repositoryRoot, "src"));

    if exist("readCbModel", "file") ~= 2
        error( ...
            "flutor:CobraDependencyMissing", ...
            "COBRA Toolbox function readCbModel was not found.");
    end
    if ~isfile(modelFile)
        error("flutor:ModelFileNotFound", "Model file not found: %s", modelFile);
    end

    model = readCbModel(modelFile);

    if strlength(boundsFile) > 0
        if ~isfile(boundsFile)
            error("flutor:BoundsFileNotFound", "Bounds file not found: %s", boundsFile);
        end
        boundChanges = readtable(boundsFile, 'TextType', 'string');
        model = flutor.applyReactionBounds(model, boundChanges);
    end

    output = flutor.runAnalysis(model, biomassReaction, struct());

    metadata = struct( ...
        "modelFile", char(modelFile), ...
        "biomassReaction", char(biomassReaction), ...
        "boundsFile", char(boundsFile), ...
        "matlabRelease", version("-release"), ...
        "matlabVersion", version, ...
        "inputSHA256", struct( ...
            "model", localSha256(modelFile), ...
            "bounds", localOptionalSha256(boundsFile)));

    [~, runName] = fileparts(modelFile);
    output.paths = flutor.writeResults( ...
        outputDirectory, ...
        runName, ...
        output, ...
        metadata);
end

function digest = localOptionalSha256(path)
    if strlength(path) == 0
        digest = "";
    else
        digest = localSha256(path);
    end
end

function digest = localSha256(path)
    fid = fopen(path, "rb");
    if fid < 0
        error("flutor:FileReadFailed", "Could not open %s.", path);
    end
    cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, "*uint8");
    clear cleanup;

    engine = java.security.MessageDigest.getInstance("SHA-256");
    engine.update(bytes);
    digest = lower(string(reshape(dec2hex(typecast(engine.digest(), "uint8"), 2).', 1, [])));
end

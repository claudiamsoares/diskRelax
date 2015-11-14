function varargout = decodeVarargin(args)
    
    MAXITER = 1e6;
    stopWhenDone = true;
    relativeimprovementsp = false;
    fixediterationsp = false;
    epsilons = [];
    savedatap = false;
    saveatiter = 0;

    for k=1:2:length(args)
        if strcmpi(args{k}, 'maxiter')
            MAXITER = args{k+1};
            stopWhenDone = false;
        elseif strcmpi(args{k},'epsilon')
            epsilons = args{k+1};
        elseif strcmpi(args{k}, 'relativeimprovements')
            epsilons = sort(args{k+1},'descend');
            relativeimprovementsp = true;
            savedatap = true;
        elseif strcmpi(args{k}, 'fixediterations')
            saveatiter = args{k+1};
            fixediterationsp = true;
            savedatap = true;
            MAXITER = max(MAXITER, saveatiter);
        end
    end

    
    if isempty(epsilons)
        epsilon = 1e-6;
    else
        epsilon = min(epsilons);
    end
    
    savedata.p = savedatap;
    savedata.relativeimprovementsp = relativeimprovementsp;
    savedata.fixediterationsp = fixediterationsp;
    savedata.epsiloncount = 1;
    savedata.itercount = 1;
    savedata.epsilonsensors = [];
    savedata.itersensors = [];
    savedata.epsiloniter = [];
    savedata.iteriter = [];
    savedata.epsilons = epsilons;
    savedata.saveatiter = saveatiter;
    
    
    varargout = {MAXITER, epsilon, stopWhenDone, savedata};
    
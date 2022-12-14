%%%%%%%%%%%%%%%%% computes the waiting times 
%%%% In: DataExp, DataPred, PosPred, tmax end of the trace for each cell, hyperparameters
%%%% Out: dt, dtc waiting times, dtc are right censored, T0 is the beginning
%%%% of the analysed region (for non-stationary signals)

function [dt,dtc,T0] = dtcomp(DataExp,DataPred,PosPred,tmax,FrameLen,EspaceInterPolyMin,Polym_speed,TaillePostMarq,TaillePreMarq,TailleSeqMarq)
    % -------Parameters----------
    FreqEchImg = (1/FrameLen); % image per second data time sampling   
    FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed); % how many interval(possible poly start position) in 1s
    Dwell = (TaillePreMarq + TailleSeqMarq + TaillePostMarq) / Polym_speed; % (s)   
    n=size(DataExp);
    nexp=n(2);
    nn=size(PosPred);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T0=[];  %%%%% will contain the start of the analyzed region
    %%%%%%%%%%%% 
    for data_i=1:nexp

    max_intensity=max(DataPred(:,data_i));
    %%%% find first hit
    ihit=min(find(DataPred(:,data_i) > max_intensity/5 ));
    if isempty(ihit)
        ihit = n(1);
    end

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    t0=(ihit-1)/FreqEchImg; %%%% t0 
    T0=[T0,t0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% compute intervals %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt=[]; dtc=[]; 
    for i=1:nn(2) %%% for all cells  
        times = find (PosPred (:,i) == 1) / FreqEchSimu - Dwell ; %%% starting times of polymerases in seconds
        %%%%%    eliminate events that don't have effect between T0 and tmax
        times = times ( times + Dwell > T0(i) & times + (TaillePreMarq)/Polym_speed < tmax(i) );  
        lt = length(times);  
        if lt > 1
            dtimes = diff(times);
            dt=[dt;dtimes];
        if tmax(i) > times(end) 
            dtc=[dtc;tmax(i)-times(end)]; %%%% the last truncated interval
        end
        else
           if lt == 1
                dtc=[dtc;tmax(i)-times(end)]; %%%% the last truncated interval
           end
           if lt==0
                dtc=[dtc;tmax(i)-T0(i) + (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed]; %%%% truncated interval
           end
        end
    end
end

numEvents = size(ripple.idx,1); %Determine number of detected spindles on the selected LFP/EEG channel
% Determine maximum y-axis values (mV) for the event filtered signals plus
% added 10%
maxY_LFPbp = max(abs(filtRipple)) + 0.1*max(abs(filtRipple));
maxY_LFP = max(abs(smooth(LFP))) + 0.1*max(abs(smooth(LFP)));

%% Calculate and plot results based on each detected ripple centered in a 1 sec window:
for i = 1:numEvents     
    % Determine how many datapoints to add to each end of the ripple to
    % make a 1 sec window:
    rippleDataPts =  ripple.idx(i,2) - ripple.idx(i,1) + 1;
	addPts = floor(1 * Fs) - rippleDataPts;
    if addPts < 0   % If the ripple is the max allowed of 300 msec duration,
        addPts = 0; %do not add any data points to the beginning and end.
    end
    start = ripple.idx(i,1) - floor(addPts/2);
    stop = ripple.idx(i,2) + ceil(addPts/2);

%     cells = size(spikeData,2);
%     spikesInSpindle = [];
    
%     for c = 1:cells                                                         % For each place cell
%         spikes = spikeData{c}/1000000;                                      % Keep all its spikes
%         spikes = spikes(spikes >= TimeStamps(start) & spikes <= TimeStamps(stop)); % Keep all its spikes that fall within the target interval
%         if ~isempty(spikes)                                                 % If there are any spikes
%             spikesInSpindle = [spikesInSpindle; [spikes c*ones(size(spikes,1),1)]]; 
%         end
%     end

%     if ~isempty(spikesInSpindle)
%         lcSpikes = lcSpikeData/1000000;
%         lcSpikes = lcSpikes(lcSpikes >= TimeStamps(start) & lcSpikes <= TimeStamps(stop));
%         lcSpikesInSpindle = [];
%         if ~isempty(lcSpikes)
%             lcSpikesInSpindle = [lcSpikesInSpindle; [lcSpikes (cells+1)*ones(size(lcSpikes,1),1)]];
%         end
%         targetIdx = rippleFiltTS >= TimeStamps(start) & rippleFiltTS <= TimeStamps(stop);
        
        ax1 = subplot(2,1,1);
        plot(ax1, timeStamps(start:stop), filtRipple(start:stop), 'c')
        hold on
        plot(ax1, timeStamps(ripple.idx(i,1):ripple.idx(i,2)),...
            filtRipple(ripple.idx(i,1):ripple.idx(i,2)), 'b')
        plot(ax1, timeStamps(start:stop), LFP(start:stop), 'Color', [0.5 0.5 0.5])
        plot(ax1, timeStamps(start:stop), rmsSignal(start:stop), 'k')
        plot(ax1, timeStamps(start:stop), upperThreshold* ones(stop-start+1,1), 'g')
        plot(ax1, timeStamps(start:stop), lowerThreshold* ones(stop-start+1,1), 'r')
        hold off
        switch ripple.state(i)
            case 2
                stage = 'SWS';
            case 3
                stage = 'REM';
            case 4
                stage = 'QW';
            case 6
                stage = 'TR ';
        end
        title(['Ripple ' num2str(i) ' of ' num2str(numEvents) '          Stage: ' stage])
        xlim([timeStamps(start) timeStamps(stop)])
        ylim([-maxY_LFPbp maxY_LFPbp])
        ylabel(ax1, 'HC LFP')
        
%         ax2 = subplot(3,1,2);
%         map=lines(7);
%         colorUnit = map(mod(spikesInSpindle(1,2),7)+1,:);
%         plot(ax2, [spikesInSpindle(1,1) spikesInSpindle(1,1)], [spikesInSpindle(1,2)-.5 spikesInSpindle(1,2)+.5], 'Color', colorUnit)
%         hold on
%         for m = 2:size(spikesInSpindle,1)
%             colorUnit = map(mod(spikesInSpindle(m,2),7)+1,:);             %de2bi(mod(spikesInSpindle(m,2),7),3);
%             plot(ax2, [spikesInSpindle(m,1) spikesInSpindle(m,1)], [spikesInSpindle(m,2)-.5 spikesInSpindle(m,2)+.5], 'Color', colorUnit)
%         end
% %         scatter(ax2,spikesInSpindle(:,1),spikesInSpindle(:,2), 'b', 'd')
% %         hold on
%         if ~isempty(lcSpikesInSpindle)
%             scatter(ax2,lcSpikesInSpindle(:,1),lcSpikesInSpindle(:,2), 'r', 's')
%         end
%         hold off
%         ylim([0 cells+1])
%         xlim([TimeStamps(start) TimeStamps(stop)])
%         ylabel(ax2, 'Place Cells')
% %         ax2.TickDirMode = 'manual';
% %         ax2.TickDir = 'out';
%         
        ax2 = subplot(2,1,2);
        plot(ax2, timeStamps(start:stop), LFP(start:stop), 'Color', [0.5 0.5 0.5])
        hold on
        plot(ax2, timeStamps(ripple.idx(i,1):ripple.idx(i,2)),...
            LFP(ripple.idx(i,1):ripple.idx(i,2)), 'b')
        hold off
        xlim([timeStamps(start) timeStamps(stop)])
        ylim([-maxY_LFP maxY_LFP]);
        ylabel(ax2, 'Hippocampal LFP')
        xlabel(ax2, 'Time(seconds)');
         pause
%     end

end
% NeuronBiomarkers
% 09/12/15
% function [biomarkers biomarkerNames] = NeuronBiomarkers()
function NeuronBiomarkers()
% Collect all neuronal biomarkers we are using
% Common format: V is a 2 column vector, V(:,1) is time, V(:,2) is voltage.

% Necessary data for biomarker calculation
% Set of traces following protocol
% Stimulus levels for each trace

% Necessary processing for biomarker calculation
% Each trace checked for AP, if AP(s) exist, indices for each AP.

% Diagnostics:

% -- Load trace file and sort into individualrdy APs --
filename = 'basicSpikes.dat';
baseTrace =  OpenNeuronTrace(filename);

assert(size(baseTrace,2) == 2);
assert(size(baseTrace,1) > 2);

% TestSplitTraceIntoAPs(baseTrace);

traces = SplitTraceIntoAPs(baseTrace);


% Sort baseTrace into individual traces

% -- Run biomarkers on each trace and compare answers to my estimates --



% --- Analyse special APs ---

% Resting AP (restingAP)
% No stimulus, for finding resting potential

% Step test
% Series of APs with different stimulus current amplitudes to find AP
% threshold, rheobase etc. Store in a big matrix.

stepIncrement = 50; % pA
stepDuration = 800; % ms
numSteps = 80;

% Ramp test
% Series of APs with different ramp currents (from 0 to stimAmp over
% rampDuration)

rampIncrement = 50; % pA
rampDuration = 500; % ms
numSteps = 80;


disp('Go!')
end

% --|| SUBFUNCTIONS ||--
function trace = OpenNeuronTrace(filename)

trace = importdata(filename);

end

% Function to split a trace into its constituent APs
function [traces numAPs] = SplitTraceIntoAPs(trace)

% Find all APs (points where trace crosses voltage threshold, excluding points within the
% time threshold caused by experimental noise)

% ///MAGIC NUMBERS\\\
threshold = 0; % Manual threshold (mV)
timeThreshold = 5; % Threshold to group crossings of threshold (ms)

% Find each instance of voltage threshold crossing from below
crossings = [];
for i = 1:length(trace(:,1))-1
    
    if trace(i,2) < threshold;
        if trace(i+1,2) > threshold;
            crossings(end+1,1) = i;
        end
    end
end

% For each crossing, remove all instances within 5 ms, leaving only the
% first crossing of the threshold
timeCrossings = trace(crossings,1);
groupedCrossings = zeros(size(crossings));
for i = 1:length(crossings)-1
    
    % If not part of a group or the first instance of crossing the threshold
    % in a group of crossings, then continue to process, otherwise move on
    if groupedCrossings(i) == 0
        nearbyCrossings = (timeCrossings(i+1:end) - timeCrossings(i)) < timeThreshold; % Future crossings occurring in time threshold
        groupedCrossings(i+1:end) = groupedCrossings(i+1:end) + nearbyCrossings; % Assign them
        assert(all(groupedCrossings < 2)); % Make sure we are properly skipping all the crossings that are counted as groups
    end
    
end
firstCrossIndices = find(groupedCrossings == 0);

% Store only the first instance from each 5 ms interval starting at first point where voltage crossed
% the threshold from below without having previously crossed in last 5 ms
% (assuming time threshold = 5)
firstCrossings = crossings(firstCrossIndices);

% Compute number of APs
numAPs = length(firstCrossings);

% If num APs is 1 or 0, return one trace, the whole trace.
% Otherwise:
if numAPs == 0
    traces = trace;
    numAPs = 0;
    return
elseif numAPs == 1
    traces = trace;
    numAPs = 1;
    return
end

% For each AP, find the minimum value of Vm before the next AP.
assert(numAPs > 1);
for i = 1:numAPs
    
    if i == 1
        startIdx(i) = 1;
    else
        startIdx(i) = endIdx(i-1)+1;
    end
    
    if i == numAPs
        endIdx(i) = length(trace);
    else
        % Calculate end of this trace
        
        % Get voltage from this AP's crossing to next AP's crossing
        curAPUpIndex = firstCrossings(i);
        nextAPUpIndex = firstCrossings(i+1);
        voltageDuringCurrentAP = trace(curAPUpIndex:nextAPUpIndex,2);
        
        % Find the minimum voltage and use it to calculate end index
        [minVm minVmIdx] = min(voltageDuringCurrentAP);
        endIdx(i) = minVmIdx + curAPUpIndex-1; % Add back in previous indices
    end
    
    % Add the AP to the cell of APs
    traces{i} = trace(startIdx(i):endIdx(i),:);
    
end

assert(length(traces) == numAPs);
end

% --Gradient--
function dVdt = VoltageGradient(trace)
dt = diff(trace(:,1));
dV = diff(trace(:,2));
dVdt = [trace(1:end-1,1) dV./dt]; % If dt is in ms and dV in mV, units will be V/s
end

% --Biomarkers--

% Resting Potential
function restingVm = RMP(trace)
V = trace;
restinVm = min(V(:,2));
end

% Ramp AP
function APs = RampAPs()

end
% AP Threshold

% Step rheobase

function out = StepRheobase(stepTestTraces,duration)
% Processes a step test to find the first AP

% for i = 1:length(stepTest)
% if Find an AP
% id = i
% break
% end
%
% out = stepIncrement * i;

end

% AP Peak
function [peak peakTime] = APPeak(trace)

[peak peakIdx] = max(trace(:,2));
peakTime = trace(peakIdx,1);

end

% AP Rise time
% Time from dVdt above a certain threshold, up to when AP peak was reached
function riseTime = RiseTime(trace,threshold)

% Get dVdt
dVdt = VoltageGradient(trace)
% Find max
[peak peakTime] = APPeak(trace)

% Find threshold (if it occurs)
foundThreshold = [];
for i = 1:length(dVdt)-1
    if dVdt(i,2) < threshold
        if dVdt(i+1,2) > threshold
            foundThreshold(end+1) = i;
        end
    end
    
end

numThresholds = length(foundThreshold);
if numThreshold == 1
    % Find time at threshold
    thresholdTime = trace(foundThreshold(1),1);
    % Find diffrence from peak time
    riseTime = peakTime - thresholdTime;
    assert(riseTime >= 0);
elseif numThresholds == 0
    riseTime = 'N/A'
elseif numThresholds > 1
    error('More than 1 threshold for rise time - APs possibly not properly separated')
    
end

end

% AP slope max/min
function [slopeMin slopeMax] = APSlopeMinMax(trace)

dVdt = VoltageGradient(trace);
slopeMin = min(dVdt(:,2));
slopeMax = max(dVdt(:,2));

end

% AP width full
function [fullWidth numUps numDowns] = APFullWidth(trace,threshold)

% Measure width of AP in ms at threshold (default from Davidson et al. 2014
% is 5 mV.

% Find intersections at threshold. If less than 1 up and 1 down
% intersection return "N/A".
% Otherwise, return the difference between earliest up and latest down
% intersection.

ups = [];
downs = [];
for i = 1:length(trace)-1
    
    % Find the ups (crossing theshold from below)
    
    % Find the downs (crossing threshold from above)
    
end

% Count them
numUps = length(ups);
numDowns = length(downs);

% Calculate width
if numUps < 1 || numDowns < 1 % Not enough crossings
    fullWidth = 'N/A'
elseif numUps == 1 && numDowns == 1 % Exactly one crossing of threshold each way
    fullWidth = trace(downs(1),1) - trace(ups(1),1);
    
elseif numUps > 1 || numDowns > 1 % More than one threshold crossing one or both ways
    
    % Ups and downs will be found in chronological order, so find earliest
    % crossing from below, and latest crossing from above to calculate full
    % width.
    earliestUp = ups(1);
    latestDown = downs(end);
    fullWidth = trace(latestDown,1) - trace(earliestUp,1);
end

end

function [amp tau] = FitAfterHyperpolarisation(trace,nextTrace)

% Fit an afterhyperpolarisation curve to the end of an AP using the period
% between the peak of the current trace and the peak of the next trace.

% Reconstruct the repolarising period
maxIdx = [];
[~, maxIdx(1)] = max(trace(:,2));
[~, maxIdx(2)] = max(nextTrace(:,2))

workingTrace = [trace(maxIdx(1):end,:); nextTrace(1:maxIdx(2),:)];

% AHP amplitude
[amp ampI] = min(workingTrace);

% AHP time constant
% Divide up and fit somehow. Need to discard upstroke.



end
% |||---- TESTS ----|||

    function TestSplitTraceIntoAPs(trace)
        
        [traces numAPs] = SplitTraceIntoAPs(trace);
        numAPs
        
        fig
        hold on
        for i = 1:numAPs
            plot(traces{i}(:,1),traces{i}(:,2),'color',rand(3,1))
        end
        
        disp('TestSplitTraceIntoAPs complete')
    end
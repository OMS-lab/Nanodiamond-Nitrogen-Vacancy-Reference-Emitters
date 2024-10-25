%% Photons class
%   Author  : Patrick Parkinson (patrick.parkinson@manchester.ac.uk)
%   Date    : 24/10/2017
%   Version : 0.2
%   License : CC-BY-SA (http://creativecommons.org/licenses/by-sa/4.0/)
%
%   This MATLAB class provides a handle-based method of dealing with the
%   photons produced by the iTCSPC method, which may stretch to several GB.
%
%   Changelog :
%       25/1/2018 : General cleanup and commenting (PWP)
%       12/4/2018 : Serialise and deserialise functions (for saving to
%                   file) (PWP)
%       24/01/2023: Switch cycles to marker
classdef photons < handle
    
    % Private (internal) variable storage
    properties 
        iCycle
        iChannel
        iTau
        iSync
        iMarker
        iPreallocate
        iPeriods = []
    end
    
    % Public (but set-blocked) variable
    properties (SetAccess=private)
        % Total photons
        n=0
        mode
    end
    
    properties
        scan = struct();
    end
    
    % Dependant variables
    properties (Dependent=true)
        % Count cycles of movement since start of experiment
        cycle
        % Photon arrival channel
        channel
        % Time since excitation pulse (T3 only)
        tau
        % Time since start of experiment (T2:ps, T3:sync pulses)
        sync
        % Position marker from translation stage
        marker
        % Time spent in each marker period
        periods
    end
    
    % Public methods
    methods
        % General get methods
        function c=get.tau(obj)
            % Time since excitation pulse
            if obj.mode == 3
                c = obj.iTau(1:obj.n);
            else
                warning('photons:T2:noTau','No tau information in T2 mode');
                c = [];
            end
        end
        
        function set.periods(obj,periods)
            obj.iPeriods = periods(periods>0)';
        end
        
        function p = get.periods(obj)
            p = obj.iPeriods;
        end
        
        function c=get.channel(obj)
            c = obj.iChannel(1:obj.n);
        end
        
        function c=get.sync(obj)
            c = obj.iSync(1:obj.n);
        end
        
        function c=get.cycle(obj)
            c = obj.iCycle(1:obj.n);
        end
        
        function c=get.marker(obj)
            c = obj.iMarker(1:obj.n);
        end
        
        function shrink(obj)
            % Reduce memory usage by shrinking allocated storage to data
            % size.
            obj.iChannel = obj.iChannel(1:obj.n);
            if obj.mode == 3
                obj.iTau     = obj.iTau(1:obj.n);
            end
            obj.iSync    = obj.iSync(1:obj.n);
            obj.iMarker  = obj.iMarker(1:obj.n);
            obj.iCycle   = obj.iCycle(1:obj.n);
            obj.iPreallocate = obj.n;
        end
        
        function add(obj,toadd)
            % Call to load photons into photon structure (old method)
            N = toadd.n;
            if N+obj.n > obj.iPreallocate
                error('Photons:preallocation:filled','Preallocation filled');
            end
            obj.iChannel(obj.n+1:obj.n+N)     = toadd.channel;
            obj.iCycle(obj.n+1:obj.n+N)       = toadd.cycle;
            if obj.n == 0
                misync = 0;
            else
                misync = max(obj.iSync);
            end
            if obj.mode == 3
                obj.iTau(obj.n+1:obj.n+N)     = toadd.tau;
                obj.iSync(obj.n+1:obj.n+N)    = toadd.sync+misync;
                obj.iMarker(obj.n+1:obj.n+N)  = toadd.marker;
            elseif obj.mode == 2
                obj.iSync(obj.n+1:obj.n+N)    = toadd.sync+maxmisync;
                obj.iMarker(obj.n+1:obj.n+N)  = toadd.marker;
            end
            % Need to add periods
            if isempty(obj.iPeriods)
                obj.iPeriods = toadd.periods;
            elseif length(toadd.periods) == length(obj.iPeriods)
                obj.iPeriods = obj.iPeriods + toadd.periods;
                error('photons:periods:unequalLength','Periods are different lengths');
            end
            obj.n = obj.n+N;
        end
        
        function s=size(obj)
            % Return actual size of data stored in MB
            if  obj.mode == 3
                s = 17 * obj.iPreallocate/1024^2;
            else
                s = 13 * obj.iPreallocate/1024^2;
            end
        end
        
        % Constructor
        function obj=photons(mode)
            if ~or(mode==3,mode==2)
                error('photons:mode:incorrect','Incorrect mode specified (2,3)');
            end
            obj.mode = mode;
        end
        
        % Preallocation
        function preallocate(obj,preallocate)
            % Allocate empty space to store photons - this speeds storage
            % significantly.
            if nargin ==0
                preallocate = 0.01*1024^3;
            end
            % Preallocation
            obj.iPreallocate = round(preallocate);
            obj.iChannel     = zeros(1,obj.iPreallocate,'uint8');
            if obj.mode==3
                obj.iTau     = zeros(1,obj.iPreallocate,'uint32');
            end
            obj.iSync    = zeros(1,obj.iPreallocate,'uint64');
            obj.iMarker  = zeros(1,obj.iPreallocate,'uint16');
            obj.iCycle   = zeros(1,obj.iPreallocate,'uint16');
            obj.n        = 0;
        end
        
        function show(obj)
            % Simple function to show photons per channel as a function of
            % arrival time.
            [hc1,e] = histcounts(obj.iSync(obj.iChannel==1));
            hc2     = histcounts(obj.iSync(obj.iChannel==2),e);
            if obj.mode == 2
                plot(e/1e12,hc1,e/1e12,hc2);
                xlabel('Time (s)');
            else
                plot(e,hc1,e,hc2);
                xlabel('Sync pulses');
            end
            ylabel('Counts');
            legend('Channel 1','Channel 2','location','best');
        end
        
        % To add new photons incrementally
        function consume(obj,ph)
            % During TTTR acquisition, consume/add photons into existing
            % structure. Needs to be fast.
            N = length(ph(:,1));
            if N+obj.n > obj.iPreallocate
                error('Preallocation filled');
            end
            range = obj.n+1:obj.n+N;
            obj.iChannel(range) = ph(:,1);
            obj.iTau(range)     = ph(:,2);
            obj.iSync(range)    = ph(:,3);
            obj.iCycle(range)   = ph(:,4);
            obj.iMarker(range)  = ph(:,5);
            obj.n = obj.n+N;
        end
        
        %% Save and load functions
        function serialise(obj,filename)
            % Test implementations for serialising to disc (saving)
            fn=fopen(filename,'w','a');
            fprintf(fn,'iTCSPC datafile\nPhotons:%d\nPeriods:%d\n',obj.n,numel(obj.periods));
            fnames = fieldnames(obj.scan);
            fprintf(fn,'==HEADER\n');
            for i=1:numel(fnames)
                if ischar(obj.scan.(fnames{i}))
                    fprintf(fn,'%s:%s\n',fnames{i},obj.scan.(fnames{i}));
                elseif ~isstruct(obj.scan.(fnames{i}))
                    fprintf(fn,'%s:%s\n',fnames{i},num2str(obj.scan.(fnames{i})));
                end
            end
            fprintf(fn,'==PERIODS\n');
            fprintf(fn,'%d,',obj.periods);
            % Create diffs
            dc = uint8(diff(obj.iCycle));
            ds = uint16(diff(obj.iSync));
            dm = diff(int32(obj.iMarker));
            mindifmarker=min(dm);
            dm = uint8(dm-mindifmarker);
            % Calculate required bits
            cycle_bits   = ceil(log2(double(max(dc))))+1;
            channel_bits = ceil(log2(double(max(obj.iChannel))));
            tau_bits     = ceil(log2(double(max(obj.iTau))))+1;
            sync_bits    = ceil(log2(double(max(ds))))+1;
            marker_bits  = ceil(log2(double(max(dm))))+1;
            tbits = obj.n*(cycle_bits +channel_bits+ tau_bits + sync_bits + marker_bits)/(8*1024^2);
            fprintf(1,'Writing %s MB file\n',int2str(tbits));
            fprintf(fn,'\n==PHOTONSBITS\nDiffCycle:%d\nChannel:%d\nTau:%d\nDiffSync:%d\nDiffMarker:%d\nMinDiffMarker:%d\n',cycle_bits,channel_bits,tau_bits,sync_bits,marker_bits,mindifmarker);
            % Photons out
            fprintf(fn,'\n==PHOTONS\n');
            fwrite(fn,dc,            sprintf('ubit%d',cycle_bits));
            fwrite(fn,obj.iChannel-1,sprintf('ubit%d',channel_bits));
            fwrite(fn,obj.iTau,      sprintf('ubit%d',tau_bits));
            fwrite(fn,ds,            sprintf('ubit%d',sync_bits));
            fwrite(fn,dm,            sprintf('ubit%d',marker_bits));
            % Done
            fclose(fn);
        end
    end
end

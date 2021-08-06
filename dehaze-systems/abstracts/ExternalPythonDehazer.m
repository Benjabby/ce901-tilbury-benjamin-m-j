classdef (Abstract) ExternalPythonDehazer < BaseDehazer
   
    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = false;
        PredictsT   = false;
    end
    
    properties (SetAccess = private)
        type;
        proc;
        binWrite;
        binRead;
        started = false;
        
%         logger;
%         outputLog = ""
    end
    
%     methods (Access = private)
%         function self = log(self, ~,event)
%             if(~isempty(event.Data)) 
%                 self.outputLog = append(self.outputLog, event.Data);
%             end
%         end
%     end
    
    methods
        function self = ExternalPythonDehazer(type, showWindow, postponeStart)
            self = self@BaseDehazer;
            
            if ~any(strcmp(type,["AOD";"PMHLD"]))
                error("Currently only AOD-Net and PMHLD-Net are supported");
            end

            self.type = type;
            
            self.proc = System.Diagnostics.Process;
            self.proc.StartInfo.FileName = 'python.exe';
            ipath = fileparts(mfilename('fullpath'));
            self.proc.StartInfo.Arguments = append(fullfile(ipath,"..","neural-nets","pyinterfacer.py"),' ',self.type);
            self.proc.StartInfo.CreateNoWindow = nargin==0 || ~showWindow;
            self.proc.StartInfo.UseShellExecute = false;
            self.proc.StartInfo.RedirectStandardOutput = true;
            self.proc.StartInfo.RedirectStandardInput = true;
            
%             if self.proc.StartInfo.CreateNoWindow
%                 self.proc.StartInfo.RedirectStandardError = true;
%                 self.proc.EnableRaisingEvents = true;
%                 self.logger = self.proc.addlistener('ErrorDataReceived',@self.log);
%             end
            
            if nargin<3 || ~postponeStart
                self.startProcess;
            end
        end
        
        function self = startProcess(self, force)
            if exist('force','var') && force
                self.proc.Close();
            end
            if ~self.started || self.proc.HasExited
                if self.proc.Start()
                    self.binWrite = System.IO.BinaryWriter(self.proc.StandardInput.BaseStream);
                    self.binRead = System.IO.BinaryReader(self.proc.StandardOutput.BaseStream);
                    self.started = true;
                end
            end
            
        end
        
        function self = killProcess(self, force)
            if exist('force','var') && force
                self.proc.Kill();
            else
                self.binWrite.Write(uint8(0));
                self.binWrite.Flush();
            end
        end
        
%         function self = newSequence(self, knowns)
%             newSequence@BaseDehazer(self, knowns);
% 
%         end
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img, ~)
            img = im2uint8(img);
            [h, w, ~] = size(img);
            
            byteSize = h*w*3*8;
            
            h = uint16(h);
            w = uint16(w);
            
            self.binWrite.Write(uint8(1));
            self.binWrite.Write(typecast(h,'uint8'));
            self.binWrite.Write(typecast(w,'uint8'));
            self.binWrite.Flush();
            
            img = permute(img, [3 2 1]);
            img = img(:);
            
            self.binWrite.Write(img);
        
            timeImage = self.binRead.ReadBytes(4); % First 4 bytes are the time taken
            timeImage = uint8(timeImage);
            timeImage = typecast(timeImage,'single');

            predImage = self.binRead.ReadBytes(byteSize); % Remaining bytes are the output image

            predImage = uint8(predImage);
            predImage = typecast(predImage,'double');
            predImage = reshape(predImage, [h, w, 3]);
            
            predT = [];
            predA = [];
            timeA = [];
            
            predImage = BaseDehazer.clip(predImage);
        end
        
        
        function delete(self)
            if ~isempty(self.proc), self.killProcess(true); end
        end
    end
    
end


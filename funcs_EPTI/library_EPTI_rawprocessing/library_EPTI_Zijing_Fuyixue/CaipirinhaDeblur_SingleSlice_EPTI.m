function [Kcorrected] = CaipirinhaDeblur_SingleSlice_EPTI(K, prot, evp , PhaseShiftBtwSimulSlices, SliceSep, slicenum, start_point)

% SMS encoding deblur for EPTI data with SMS multiband acquisition for one
% slice data
% This version supports 1/2 1/3 1/4 FOV shift
% By Fuyixue Wang, 2020

% Make it works for single-shot EPTI data, which includes self-navigator at
% the beginning (start point: the number of lines of self-navigaotor - 1)
% By Zijing Dong, MGH, 2021

if size(prot.sSliceArray,2) > 1;
    
    % When the value for these coordinate is very small, e.g.
    % sPosition_dSag = 1.343452e-9 then the read_meas_dat function that readin the data would not
    % recognize it and will leave the array empty so fix it here
    if isempty(prot.sSliceArray(1).sPosition_dSag) && isempty(prot.sSliceArray(2).sPosition_dSag)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dSag = 0;
            prot.sSliceArray(count).sNormal_dSag = 0;
        end
    end
    if isempty(prot.sSliceArray(1).sPosition_dCor) && isempty(prot.sSliceArray(2).sPosition_dCor)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dCor = 0;
            prot.sSliceArray(count).sNormal_dCor = 0;
        end
    end
    if isempty(prot.sSliceArray(1).sPosition_dTra) && isempty(prot.sSliceArray(2).sPosition_dTra)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dTra = 0;
            prot.sSliceArray(count).sNormal_dTra = 0;
        end
    end
    
    NormalVec = [prot.sSliceArray(1).sNormal_dSag, prot.sSliceArray(1).sNormal_dCor, prot.sSliceArray(1).sNormal_dTra].';   
    Pos(:,1) = [prot.sSliceArray(1:end).sPosition_dSag].';
    Pos(:,2) = [prot.sSliceArray(1:end).sPosition_dCor].';
    Pos(:,3) = [prot.sSliceArray(1:end).sPosition_dTra].';
    SlicePos = Pos*NormalVec;% +4.3/2;
else
    keyboard('only single slice data so cant determine correctionFac if Gz is rotated')
end
      
% SlicePos  = SlicePos(end:-1:1); 
PhaseShiftPerMM = (PhaseShiftBtwSimulSlices/SliceSep);

%Kcorrected = single(zeros(size(K)));
Kcorrected = K;
if PhaseShiftBtwSimulSlices ~= 0
    for SlcCount = 1: size(K,10)
        if abs(PhaseShiftBtwSimulSlices) == pi % FOV/2 shift
            if start_point>1
            Kcorrected(:,1:start_point-1,:,:,:,:,:,:,:,SlcCount) =  K(:,1:start_point-1,:,:,:,:,:,:,:,SlcCount)*exp(+i*PhaseShiftPerMM*(SlicePos(slicenum))/2);
            end
            Kcorrected(:,start_point:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point:2:end,:,:,:,:,:,:,:,SlcCount)*exp(+i*PhaseShiftPerMM*(SlicePos(slicenum))/2);
%             Kcorrected(:,2:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:2:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,start_point+1:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point+1:2:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*(SlicePos(slicenum))/2);
%             Kcorrected(:,2:2:end,:,:,:,:,:,:,:,SlcCount) =  K(:,2:2:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*(SlicePos(SlcCount)-pi));

        elseif abs(PhaseShiftBtwSimulSlices) == 2*pi/3 % FOV/3 shift
            Kcorrected(:,start_point:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point:3:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,start_point+1:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point+1:3:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(slicenum));
            Kcorrected(:,start_point+2:3:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point+2:3:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(slicenum));
        elseif abs(PhaseShiftBtwSimulSlices) == pi/2 % FOV/4 shift
            Kcorrected(:,start_point:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point:4:end,:,:,:,:,:,:,:,SlcCount);
            Kcorrected(:,start_point+1:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point+1:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*PhaseShiftPerMM*SlicePos(slicenum));
            Kcorrected(:,start_point+2:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point+2:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*2*PhaseShiftPerMM*SlicePos(slicenum));
            Kcorrected(:,start_point+3:4:end,:,:,:,:,:,:,:,SlcCount) =  K(:,start_point+3:4:end,:,:,:,:,:,:,:,SlcCount)*exp(-i*3*PhaseShiftPerMM*SlicePos(slicenum));
        end
    end
end




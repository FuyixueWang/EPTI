%**************************************************************************%
function [Rdc, Pdc, Sdc] = mrir_measdat_vox2ras__dircos(prot)
%MRIR_MEASDAT_VOX2RAS__DIRCOS
%
% varargout = mrir_measdat_vox2ras(varargin)

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/17
% $Id: mrir_measdat_vox2ras.m,v 1.2 2014/03/03 16:00:56 jonp Exp $
%**************************************************************************%

  Sdc = [prot.sSliceArray(1).sNormal_dSag;
         prot.sSliceArray(1).sNormal_dCor;
         prot.sSliceArray(1).sNormal_dTra];

  inplanerot = prot.sSliceArray(1).dInPlaneRot;


  %==--------------------------------------------------------------------==%

  orientation = mrir_measdat_vox2ras__principalorientation(Sdc(3),Sdc(1),Sdc(2));

  % in-plane rotation is rotation about the slice-normal direction, so
  % ignoring that for now, there can be rotation either about the PE
  % direction (which by definition does not change the PE vector) or about
  % the RO direction, which we consider first. since PE is orthogonal to RO,
  % the PE direction can be found by

  switch(orientation),

   case 'SMS_TRANSVERSE', % PE dir is COR (A>>P), RO dir is SAG
    Pdc = cross(Sdc/norm(Sdc), [+1;0;0]);
   case 'SMS_CORONAL',    % PE dir is SAG (R>>L), RO dir is TRA
    Pdc = cross(Sdc/norm(Sdc), [0;0;+1]);
   case 'SMS_SAGITTAL',   % PE dir is COR (A>>P), RO dir is TRA
    Pdc = cross(Sdc/norm(Sdc), [0;0;-1]);
  end;

  % direction cosines assumed to be unit normal, and ||axb|| = ||a|| ||b|| sin(theta)
  Pdc = Pdc/norm(Pdc);

  % calculate the readout direction as GR = GS x GP
  Rdc = cross(Sdc, Pdc);

  % direction cosines assumed to be unit normal
  Rdc = Rdc/norm(Rdc);

  % apply the inplane rotation, if its not zero
  if ( inplanerot ~= 0.0 ),
    Pdc(1) = cos(inplanerot)*Pdc(1) - sin(inplanerot)*Rdc(1);
    Pdc(2) = cos(inplanerot)*Pdc(2) - sin(inplanerot)*Rdc(2);
    Pdc(3) = cos(inplanerot)*Pdc(3) - sin(inplanerot)*Rdc(3);

    % recalculate the readout direction
    Rdc = cross(Sdc, Pdc);

  end;

  return;


%**************************************************************************%



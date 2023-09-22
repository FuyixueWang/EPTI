function rval = mrir_measdat_vox2ras__principalorientation(ntra, nsag, ncor)
%MRIR_MEASDAT_VOX2RAS__PRINCIPALORIENTATION
%
% varargout = mrir_measdat_vox2ras(varargin)
%
%
% see also MRIR_MEASDAT_VOX2RAS_ALT.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/17
% $Id: mrir_measdat_vox2ras.m,v 1.2 2014/03/03 16:00:56 jonp Exp $
%**************************************************************************%

  abs_sag = abs(nsag);
  abs_cor = abs(ncor);
  abs_tra = abs(ntra);

  almost_equal_sc = isAlmostEqual(abs_sag, abs_cor);
  almost_equal_st = isAlmostEqual(abs_sag, abs_tra);
  almost_equal_tc = isAlmostEqual(abs_tra, abs_cor);

  % mainly transverse case, which captures most special cases
  if (   (almost_equal_sc && almost_equal_st) ...
         || (almost_equal_sc && (abs_sag < abs_tra)) ...
         || (almost_equal_st && (abs_sag > abs_cor)) ...
         || (almost_equal_tc && (abs_cor > abs_sag)) ...
         || ((abs_sag > abs_cor) && (abs_sag < abs_tra)) ...
         || ((abs_sag < abs_cor) && (abs_cor < abs_tra)) ...
         || ((abs_sag < abs_tra) && (abs_tra > abs_cor)) ...
         || ((abs_cor < abs_tra) && (abs_tra > abs_sag))),
    rval = 'SMS_TRANSVERSE';
  elseif ( (almost_equal_sc && (abs_sag > abs_tra)) ...
           || (almost_equal_st && (abs_sag < abs_cor)) ...
           || ((abs_sag < abs_cor) && (abs_cor > abs_tra)) ...
           || ((abs_sag > abs_tra) && (abs_sag < abs_cor)) ...
           || ((abs_sag < abs_tra) && (abs_tra < abs_cor)))
    rval = 'SMS_CORONAL';
  elseif ( (almost_equal_tc && (abs_cor < abs_sag)) ...
           || ((abs_sag > abs_cor) && (abs_sag > abs_tra)) ...
           || ((abs_cor > abs_tra) && (abs_cor < abs_sag)) ...
           || ((abs_cor < abs_tra) && (abs_tra < abs_sag)))
    rval = 'SMS_SAGITTAL';
  else
    error('invalid slice orientation encountered');
  end;

  return;


%**************************************************************************%
function bool = isAlmostEqual(val_a, val_b)
%ISALMOSTEQUAL
%
% varargout = isAlmostEqual(varargin)
%
%
% see also ISALMOSTEQUAL_ALT.

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/nov/17
% $Id: mrir_measdat_vox2ras.m,v 1.2 2014/03/03 16:00:56 jonp Exp $
%**************************************************************************%

  diff = val_a - val_b;

  if ( (diff >= -1.0e-6) && (diff <= 1.0e-6) ),
    bool = 1;
  else,
    bool = 0;
  end;

  return;



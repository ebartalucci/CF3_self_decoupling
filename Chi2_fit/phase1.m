%############################################################################
%#
%#                          function PHASE
%#
%#      applies phase-correction to complex, 1D data
%#
%#      usage: pdata=phase(data,phase);
%#         or: pdata=phase(data,[-1.23,0.001,123]);
%#
%#             phase = values for phase-correction (see phasetool.m for details)
%#     
%#      (c) P. Blümler 1/03
%############################################################################
%----------------------------------------------------------------------------
%  version 1.1 PB 21/1/03    (please change this when code is altered)
%----------------------------------------------------------------------------


  function result=phase(data,pval)
  data=squeeze(data);
  dim=dimension(data);
  
  if dim ~= 1
      errordlg('ERROR: data input array is NOT ONEDIMENSIONAL!');
      return
  end
  if nargin ~= 2
      errordlg('ERROR: not all input data specified!');
      return
  end
  ns=size(data);
  if ns(2)==1
      data=data';
  end
  n_points=length(data);
  x=linspace(1,n_points,n_points)-pval(3);
  phi=pval(1)/180*pi+pval(2)*x/1000;
  phas_vector=complex(cos(phi),-sin(phi));
  result=data.*phas_vector;

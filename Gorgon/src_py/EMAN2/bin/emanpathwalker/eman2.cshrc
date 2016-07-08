#!/bin/csh
setenv EMAN2DIR /Applications/EMAN2
set path=(${EMAN2DIR}/bin ${EMAN2DIR}/extlib/bin $path)
if ( $?PYTHONPATH ) then
else
setenv PYTHONPATH
endif
setenv PYTHONPATH ${EMAN2DIR}/lib:${EMAN2DIR}/bin:${EMAN2DIR}/extlib/site-packages:${EMAN2DIR}/extlib/site-packages/ipython-1.2.1-py2.7.egg:${PYTHONPATH}

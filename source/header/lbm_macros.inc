! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Konstantin Kleinheinz <k.kleinheinz@grs-sim.de>
! Copyright (c) 2011-2012 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012, 2014-2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2013-2015, 2019 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF SIEGEN “AS IS” AND ANY EXPRESS
! OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
! OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL UNIVERSITY OF SIEGEN OR CONTRIBUTORS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
??! We include two accesses to the pdf variables
??! - array of structures (AOS)
??! - structure of arrays (SOA)
??!
??IFNDEF (SOA) THEN
?? LOGICAL :: SOA = .false.
??ENDIF
??!
??! Streaming layout
??!   PUSH:  read element-local pdf and store to adjacent element
??!   PULL:  read adjacent element pdf and store to local element
??! Multi-level only works with Pull
??!
??!
??IFNDEF (PUSH) THEN
?? LOGICAL :: PUSH = .false.
??ENDIF
??! Macro to get neighbor element index from current element index, connectivity
??! array and direction
?? macro :: NghElemIDX(neigh, iDir, stencil, index, QQ, nElems) &
??   & = ?ElemIDX?(?neigh?(?NGPOS?(?NgDir?(?iDir?, ?stencil?), ?index?, ?nElems?)), &
??   & ?QQ?, ?nElems?)
??!
??!
??! NgDir: The ELEMENT direction with which the streaming operation takes places
??!        Mainly used in the reconstruct_connectivity routine to build the neighbor array.
??!        Also used in boundary routine
??! nScalars: number of scalars in state array
??! QQ: number of stencil directions
?? if(PUSH) then
?? macro :: FETCH(iDir,iField,node,QQ,nScalars,nElems)       &
??   &        = ?IDX?( ?iDir?+(?iField?-1)*?QQ?, ?node?, ?nScalars?, ?nElems?)
?? macro :: SAVE( iDir,iField,node,QQ,nScalars,nElems,neigh) &
??   &        = ?neigh?( ?NGPOS?(?iDir?,?node?,?nElems?) ) + ?NGOFFSET?(?iField?, ?QQ?, ?nElems?) + ?nScalars?*0
?? macro :: StreamName = push
?? macro :: NgDir(iDir) = ?iDir?
?? macro :: InvNgDir(iDir, stencil) = ?stencil?%cxDirInv( ?iDir? )
?? macro :: post2pre( omega ) = ( 1.0_rk + ?omega?*0.0_rk )
?? macro :: pre2post( omega ) = ( 1.0_rk + ?omega?*0.0_rk )
?? else !PULL
?? macro :: FETCH(iDir,iField,node,QQ,nScalars,nElems,neigh) &
??   &        = ?neigh?(?NGPOS?(?iDir?,?node?,?nElems?))+?NGOFFSET?(?iField?,?QQ?,?nElems?)+?nScalars?*0
?? macro :: SAVE( iDir,iField,node,QQ,nScalars,nElems)      &
??   &        = ?IDX?( ?iDir?+(?iField?-1)*?QQ?, ?node?, ?nScalars?, ?nElems?)
?? macro :: StreamName = pull
?? macro :: NgDir(iDir, stencil) = ?stencil?%cxDirInv( ?iDir? )
?? macro :: InvNgDir(iDir) = ?iDir?
?? macro :: post2pre( omega ) = ( 1.0_rk - ?omega? )
?? macro :: pre2post( omega ) = ( 1.0_rk - ?omega? )
?? endif
??
??! Access to the Connectivity Array
?? macro :: NGPOS(iDir,node,nElems) = (?iDir?-1)*?nElems? + ?node?
??!
?? IF(SOA) THEN
??!--------------------------------------------
??!    S O A - Structure of arrays layout
??!-------------------------------------------
??! These are the array index access macros
??! Access to the state array
?? macro :: IDX(varPosiDir,node, nScalars, nElems) = (?varPosiDir?-1)*?nElems? + ?node? + ?nScalars?*0
??! Offset in NGPOS
?? macro :: NGOFFSET( iField, QQ, nElems ) = (?iField?-1)*?QQ?*?nElems?
??! Position in Element list
?? macro :: ElemIDX(index, QQ, nElems) = mod(?index?-1,?nElems?) + 1 + ?QQ?*0
??! Direction
?? macro :: DirIDX(index, QQ, nElems)  = int((?index?-1)/?nElems?) + 1 + ?QQ?*0
??!
?? ELSE
??!--------------------------------------------
??!    A O S - Array of structures layout
??!-------------------------------------------
??! These are the Array index access macros
??! Access to the State Array
?? macro :: IDX(varPosiDir,node,nScalars) = (?node?-1)*?nScalars? + ?varPosiDir?
??! Offset after NGPOS
?? macro :: NGOFFSET( iField, QQ ) = (?iField?-1)*?QQ?
??! Position in Element list
?? macro :: ElemIDX(index, QQ) = int((?index?-1)/?QQ?) + 1
??! Direction
?? macro :: DirIDX(index, QQ, nElems)  = mod(?index?-1,?QQ?) + 1 + ?nElems?*0
?? ENDIF
??!
?? text :: dir_vector
!CDIR NODEP
!IBM* INDEPENDENT
!DIR$ IVDEP
!DIR$ VECTOR aligned NONTEMPORAL(outState)
?? end text dir_vector
??!
?? text :: dir_novec
!CDIR NODEP
!IBM* NOVECTOR
!DIR$ NOVECTOR
?? end text dir_novec

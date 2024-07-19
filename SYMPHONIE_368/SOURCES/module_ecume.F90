      module module_ecume
!______________________________________________________________________
! SYMPHONIE ocean model
! release 362 - last update: 06-01-23
!______________________________________________________________________
!...............................................................................
! Version Date      Description des modifications
! v362    06-01-23  mise en service
!...............................................................................
!    _________                    .__                  .__             !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................
      INTEGER &
        JITER &
       ,NITERFL &
       ,JCV   &
       ,NITERMAX &
       ,NITERSUP &
       ,JJ


      LOGICAL &
        OPCVFLX &
       ,LPRECIP &
       ,LPWEBB

      REAL &
        PSST  &
       ,PPS   &
       ,PPA   &
       ,ZVMOD &
       ,PTA   &
       ,PQA   &
       ,ZSSS  &
       ,PEXNA &
       ,PEXNX &
       ,ZFOES &
       ,PQSAT &
       ,ZPA   &
       ,ZFOESA &
       ,ZQSATA &
       ,XG    &
       ,XP00  &
       ,XKARMAN &
       ,PZREF &
       ,PUREF &
       ,XZ0 &
       ,XTT  &
       ,ZUTU &
       ,ZUTT &
       ,ZUTQ &
       ,ZETV &
       ,ZLMOMIN &
       ,ZLMOMAX &
       ,ZLMOU &
       ,ZLMOT &
       ,ZBTA &
       ,ZGMA &
       ,XPI &
       ,ZSQR3 &
       ,Z0TSEA &
       ,Z0QSEA &
       ,PZ0SEA &
       ,PZ0HSEA &
       ,PCHN &
       ,PCEN &
       ,ZRDSRV &
       ,ZLOGUS10,ZLOGTS10 &
       ,ZCHIK &
       ,ZCHIC &
       ,ZPSIK &
       ,ZPSIC &
       ,ZPSI_U &
       ,ZPSI_T &
       ,ZPSIU &
       ,ZPSIT &
       ,ZDU  &
       ,ZDT  &
       ,ZDQ  &
       ,ZLVA &
       ,ZLVS &
       ,ZCPA &
       ,ZVISA &
       ,ZDDU &
       ,ZDDT &
       ,ZDDQ &
       ,ZUSR &
       ,ZTSR &
       ,ZQSR &
       ,ZUSR0 &
       ,ZTSR0 &
       ,ZQSR0 &
       ,ZDUSTO &
       ,ZDTSTO &
       ,ZDQSTO &
       ,ZDELTAU10N &
       ,ZDELTAT10N &
       ,ZDELTAQ10N &
       ,ZPARUN,ZPARTN,ZPARQN & 
       ,ZORDOU,ZORDOT,ZORDOQ &
       ,ZCDIRU,ZCDIRT,ZCDIRQ &
       ,ZDUSR0,ZDTSR0,ZDQSR0 &
       ,ZP00 &
       ,ZWW  &
       ,ZEFWEBB &
       ,XUNDEF &
       ,PRHOA  &
       ,PCE  &
       ,PEXNS &
       ,PRESA &
       ,PVMOD &
       ,PCD,PCH,PCDN,PRI &
       ,PUSTAR,PSFTH,PSFTQ,PRAIN,ZTAU,ZHF,ZEF,ZTAUR,ZRF,ZAC &
       ,ZTAC,ZCPWA,ZDQSDT,ZDWAT,ZDTMP,ZBULB,ZUSTAR2,PTAU


       REAL,DIMENSION(0:5):: &
        ZCOEFU &
       ,ZCOEFT &
       ,ZCOEFQ 

      end module module_ecume

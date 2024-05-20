subroutine write_record(elat,elon,rlat,rlon,stn,ntw,dt,nps,seismogram)
      !
            implicit none
      !
            include "spherical.h"
      !
            integer nps
            real elat,elon,rlat,rlon
            double precision dt,seismogram(NCOMP,NTMAX)
            character*8 stn
            character*2 ntw
      !
            integer i,loutdir
            integer lnblnk
            real ar,at,ap,ag, agrv, ahrr, ahrt, ahrp, ahtt, ahpp, ahtp
            double precision rot1,rot2
            character*12 chn
            character*80 outdir
            character*80 newfile,recordfilez,recordfilee,recordfilen
      !
      !       parameters for rotating from longitudinal and transverse components
      !       to N-S and E-W components
      !
            call rotate(elat,elon,rlat,rlon,rot1,rot2)
      !
            outdir='output'
            loutdir=lnblnk(outdir)
      !
            chn='LHZ'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilez=outdir(1:loutdir)//'/'//newfile
            open(unit=13,file=recordfilez)
      !
            chn='LHN'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilen=outdir(1:loutdir)//'/'//newfile
            open(unit=14,file=recordfilen)
      !
            chn='LHE'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=15,file=recordfilee)
      !
            chn='GRV'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=16,file=recordfilee)
      
            chn='LHZC.PGRAV'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=17,file=recordfilee)
      
            ! Vertical strain component 
            chn='LHZZ.STRAIN'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=18,file=recordfilee)

            ! Radial component (as in radial and transverse, not R THETA PHI)
            chn='LHRR.STRAIN'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=19,file=recordfilee)
      
            ! Transverse component
            chn='LHTT.STRAIN'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=20,file=recordfilee)
      
            ! Vertical/Radial component  (as in radial and transverse, not R THETA PHI)
            chn='LHZR.STRAIN'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=21,file=recordfilee)

            ! Vertical/transverse (as in radial and transverse, not R THETA PHI)
            chn='LHZT.STRAIN'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=22,file=recordfilee)
      
            ! Radial/transverse component
            chn='LHRT.STRAIN'
            call mnam(stn,ntw,chn,'nmsyn',outdir,loutdir,newfile)
            recordfilee=outdir(1:loutdir)//'/'//newfile
            open(unit=23,file=recordfilee)
      


      
      !
            do i=1,nps
              ar=sngl(seismogram(1,i))
              at=-sngl(seismogram(2,i)*rot1-seismogram(3,i)*rot2)
              ap=sngl(seismogram(2,i)*rot2+seismogram(3,i)*rot1)
              ag=sngl(seismogram(4,i))
              agrv=sngl(seismogram(5,i))
              ! gravity strain components
              ahrr=sngl(seismogram(6,i))
              ahtt=sngl(seismogram(7,i))
              ahpp=sngl(seismogram(8,i))
              ahrt=sngl(seismogram(9,i))
              ahrp=sngl(seismogram(10,i))
              ahtp=sngl(seismogram(11,i))
      
              write(13,*) sngl(dble(i-1)*dt),ar
              write(14,*) sngl(dble(i-1)*dt),at
              write(15,*) sngl(dble(i-1)*dt),ap
              write(16,*) sngl(dble(i-1)*dt),ag
      !        write(13,*) sngl(dble(i-1)*dt),sngl(seismogram(1,i))
      !        write(14,*) sngl(dble(i-1)*dt),sngl(seismogram(2,i))
      !        write(15,*) sngl(dble(i-1)*dt),sngl(seismogram(3,i))
      !        write(16,*) sngl(dble(i-1)*dt),sngl(seismogram(4,i))
              ! Grav acc 
              write(17,*) sngl(dble(i-1)*dt),agrv
              ! Grav strain 
              write(18,*) sngl(dble(i-1)*dt),ahrr 
              write(19,*) sngl(dble(i-1)*dt),ahtt 
              write(20,*) sngl(dble(i-1)*dt),ahpp 
              write(21,*) sngl(dble(i-1)*dt),ahrt 
              write(22,*) sngl(dble(i-1)*dt),ahrp 
              write(23,*) sngl(dble(i-1)*dt),ahtp 

            enddo
      !
            close(unit=13)
            close(unit=14)
            close(unit=15)
            close(unit=16)
            close(unit=17)
            close(unit=18)
            close(unit=19)
            close(unit=20)
            close(unit=21)
            close(unit=22)
            close(unit=23)
      !
            return
            end
      !===============================================================================
      
            subroutine rotate(elat,elon,rlat,rlon,rot1,rot2)
      !
      !       rot1 and rot2 rotate longitudinal and transverse components
      !       to N-S and E_W components (based upon Dahlen & Tromp eqn 10.3)
      !
            implicit none
      !
            include "spherical.h"
      !
            real elat,elon,rlat,rlon
            double precision rot1,rot2
      !
            double precision thetar,thetas,phir,phis,theta,sint
            double precision geoco
            data geoco/0.993277d0/
      !
      !       convert to a geocentric colatitude and longitude in radians
      !
            thetas=PI/2.0d0-datan(geoco*dtan(dble(elat)*PI/180.0d0))
            thetar=PI/2.0d0-datan(geoco*dtan(dble(rlat)*PI/180.0d0))
            phis=dble(elon)*PI/180.0d0
            phir=dble(rlon)*PI/180.0d0
      !
            theta=dacos(dcos(thetar)*dcos(thetas)+dsin(thetar)*dsin(thetas)*dcos(phir-phis))
            sint=dsin(theta)
            if(sint.lt.1.0d-09) sint=1.0d-09
            rot1=(dsin(thetar)*dcos(thetas)-dcos(thetar)*dsin(thetas)*dcos(phir-phis))/sint
            rot2=dsin(thetas)*dsin(phir-phis)/sint
      !
            return
            end
      !===============================================================================
      
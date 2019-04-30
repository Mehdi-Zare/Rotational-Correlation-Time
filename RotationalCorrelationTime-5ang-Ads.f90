program fun
implicit none

integer  ,  parameter                           ::      dp = selected_real_kind(15, 307)
character(len=50)                               ::      file1 
character(len=50)                               ::      file2 
character(len=50)                               ::      filename
integer                                         ::      error_flag, alloc_err                   
integer                                         ::      iii,jjj,kkk,nline1, nline2, totlines
integer,   parameter                            ::      metalMM=1229
integer,   parameter                            ::      metalQM=51
integer,   parameter                            ::      Ads=10
integer,   parameter                            ::      Water=2200*3
real(dp), dimension(2)                          ::      Xox, Yox, Zox  ! Oxygen atom coordinates            
real(dp), dimension(6600)                       ::      Xwi,Ywi,Zwi,Xwf,Ywf,Zwf, dipoli, dipolf        ! i is for initial and f for final
real(dp), dimension(10)                         ::      Xai,Yai,Zai,Xaf,Yaf,Zaf
character(len=256), dimension(6600)             ::      waterID
character(len=256), dimension(10)               ::      AdsID
real(dp)                                        ::      a               ! for FLAG
real(dp)                                        ::      timefunction
real(dp),allocatable, dimension(:,:)            ::      RESULTS
integer                                         ::      sizeofRESULTS, num
real(dp)                                        ::      ItoR, mean, SD, sumval
real(dp)                                        ::      timeint, maxtime
! Get the number of Images we want from the user
!write(*,*) "PLEASE insert the number of structures you have, they need to start with image-00001"
!read(*,*) sizeofRESULTS
!write(*,*) " WHAT IS TIME INTERVAL BETWEEN YOUR IMAGES? Insert IN pico-second   "
!read(*,*) timeint
!write(*,*) " WHAT is the maximum time you want to have Cmio at ? Insert in pico-second "
!read(*,*) maxtime
sizeofRESULTS=100
timeint=0.05_dp
maxtime=1.0_dp

sizeofRESULTS=sizeofRESULTS-1
allocate(RESULTS(INT(maxtime/timeint)+1,sizeofRESULTS))
RESULTS=0.0_dp

write(*,*) " The program STARTS" 
do jjj=1,sizeofRESULTS
  write(file1,100) jjj
  100 format('image-', I5.5)
 write(*,*) " The first image is " , jjj
  kkk=1  ! controls writing data in RESULTS array from 1 to the end
  do iii=jjj,jjj+INT(maxtime/timeint)
    if (iii > sizeofRESULTS) exit               
    write(file2,200) iii+1
    200 format('image-', I5.5)
!    write(*,*) file2
    call countline(file1,nline1)          
         if ( error_flag == 1  ) then
            call EXIT(0)           !FLAG
         end if
    call countline(file2,nline2)
         if ( error_flag == 1  ) then
            call EXIT(0)           !FLAG
         end if
    totlines=metalMM+Water+metalQM+Ads+2
         if ( totlines /= nline1 ) then
                  write(*,*) "Total lines in ", file1, " must be  ", totlines
                  call EXIT(0)    !FLAG
          else if ( totlines /= nline2 ) then
                  write(*,*) "Total lines in ", file2, " must be  ", totlines
                  call EXIT(0)    !FLAG
          end if 
  call ReadConfig(file1,metalMM,metalQM,Ads,Water,waterID,Xwi,Ywi,Zwi,AdsID,Xai,Yai,Zai)
  call ReadConfig(file2,metalMM,metalQM,Ads,Water,waterID,Xwf,Ywf,Zwf,AdsID,Xaf,Yaf,Zaf)
  call RotationalCorrelationFunction(Ads,Water,waterID,Xwi,Ywi,Zwi,AdsID,Xai,Yai,Zai,Xwf,Ywf,Zwf,Xaf,Yaf,Zaf,timefunction)

  RESULTS(kkk,jjj)=timefunction
  kkk=kkk+1
  end do
end do

ItoR=REAL(sizeofRESULTS)
!write(*,*) "real is " , ItoR
!collect data and avarage them over different time steps
!write(*,*) RESULTS(1,:)   
open    ( unit = 5000, file = 'RESULTS', status = 'new')
write   ( 5000, 5010)
5010 format ( 3x, "t(ps)", 11x, "Cmio(t)", 16x,"data-points", 16x, "St.dev")
write   (5000,5020)
5020 format (2x, "=============================================",&
                "================================================")
do iii=1,INT(maxtime/timeint)+1
  ! Calculate Mean
  num=0
  sumval=0.0_dp
  do jjj=1,sizeofRESULTS
        if (RESULTS(iii,jjj) /= 0.0_dp) then
          sumval=sumval+RESULTS(iii,jjj)
          num=num+1
        end if
  end do
  ItoR=REAL(num)
 ! write(*,*) "SUM IS ", sumval, "OR IS ", SUM(RESULTS(iii,:))
  mean=sumval/ItoR
  !write(*,*) "MEAN is " , mean

  ! Calculate Standard deviation
  sumval=0.0_dp
  do jjj=1,sizeofRESULTS
        if (RESULTS(iii,jjj) /= 0.0_dp) then
          write(*,*) "VALUE IS " , RESULTS(iii,jjj)  !FLAG
          sumval=sumval+(RESULTS(iii,jjj)-mean)**2
        end if
  end do
  SD=SQRT((1.0_dp/(num-1))*(sumval))
  !write(*,*) "NUM=" ,num
  !write(*,*) "SD=" , SD
write   (5000, 5030) iii*timeint, mean, num, SD
5030 format (2x, F6.2, 3x , ES20.13, 12x, I5, 13x, ES20.13)

end do
close(5000)


!write(*,*) size(RESULTS)
!a=Xaf(10)-Xai(10)   !FLAG
!write(*,*) a        !FLAG
deallocate(RESULTS,stat = alloc_err)

end program fun


subroutine RotationalCorrelationFunction(Ads,Water,waterID,Xwi,Ywi,Zwi,AdsID,Xai,Yai,Zai,Xwf,Ywf,Zwf,Xaf,Yaf,Zaf,timefunction)
implicit none

integer  ,  parameter                           ::      dp = selected_real_kind(15, 307)

!local variable
integer                                         ::      ierror,  error_flag, alloc_err
integer                                         ::      ii, jj, kk
real(dp), dimension(2)                          ::      Xox, Yox, Zox  ! Oxygenatom coordinates
real(dp)                                        ::      distance1i, distance2i, distance1f, distance2f, unitdipol
real(dp)                                        ::      dotdipolXi, dotdipolYi, dotdipolZi, dotdipolXf
real(dp)                                        ::      dotdipolYf,dotdipolZf,lengthi,lengthf
real(dp)                                        ::      Xi,Yi,Zi, Xf, Yf, Zf
real(dp), dimension(Water)                      ::      dipoli, dipolf
!integer,   parameter                            ::      metalMM=1229
!integer,   parameter                            ::      metalQM=51
!integer,   parameter                            ::      Ads=10
!integer,   parameter                            ::      Water=2200

!parameter types and definition
integer,                              intent(in)                  :: Ads,Water
real(dp), dimension(Water),           intent(in)                  :: Xwi,Ywi,Zwi,Xwf,Ywf,Zwf
character(len=256),dimension(Water),  intent(in)                  :: waterID
real(dp), dimension(Ads),             intent(in)                  :: Xai,Yai,Zai,Xaf,Yaf,Zaf
character(len=256),dimension(Ads),    intent(in)                  :: AdsID
real(dp),                             intent(out)                 :: timefunction

! Get the Oxygen atom coordinate of the adsorbate, it is the same for all
! conformations

jj=1 
do ii=1,Ads
        if (AdsID(ii) == 'O' ) then
        Xox(jj)=Xai(ii); Yox(jj)=Yai(ii); Zox(jj)=Zai(ii)
        jj=jj+1
        end if
end do

!write(*,*) Xox, Yox, Zox  !FLAG

dipoli=0.0_dp
dipolf=0.0_dp
kk=1 ! controls the number of water molecule in 5 ang at time zero (initial)
jj=1 ! control the number of water molecules that still exist in 5 ang at time t (final)
do ii=1,Water
        if (waterID(ii) == 'Ow') then
                distance1i=SQRT(((Xwi(ii)-Xox(1))**2) + ((Ywi(ii)-Yox(1))**2)  + ((Zwi(ii)-Zox(1))**2) )
                distance2i=SQRT(((Xwi(ii)-Xox(2))**2) + ((Ywi(ii)-Yox(2))**2)  + ((Zwi(ii)-Zox(2))**2) )
                if (distance1i <= 5 .or. distance2i <= 5) then
 !                       write(*,*) "distances are ", distance1i, distance2i
                        dotdipolXi=(((Xwi(ii+1)+Xwi(ii+2))/2.0_dp)-Xwi(ii))
                        dotdipolYi=(((Ywi(ii+1)+Ywi(ii+2))/2.0_dp)-Ywi(ii))
                        dotdipolZi=(((Zwi(ii+1)+Zwi(ii+2))/2.0_dp)-Zwi(ii))
                        lengthi=SQRT(dotdipolXi**2+dotdipolYi**2+dotdipolZi**2)
                        Xi=((dotdipolXi/lengthi)**2)
                        Yi=((dotdipolYi/lengthi)**2)
                        Zi=((dotdipolZi/lengthi)**2)
                        unitdipol=Xi+Yi+Zi
                        dipoli(kk)=unitdipol
!                       write(*,*) "ii= ",ii, 'dopoli ', dipoli(kk)
                        kk=kk+1
!                       write(*,*) "ii= " , ii, "and distance1i is" , distance1i, "distance2i= " , distance2i
                        distance1f=SQRT(((Xwf(ii)-Xox(1))**2) + ((Ywf(ii)-Yox(1))**2)  + ((Zwf(ii)-Zox(1))**2) )
                        distance2f=SQRT(((Xwf(ii)-Xox(2))**2) + ((Ywf(ii)-Yox(2))**2)  + ((Zwf(ii)-Zox(2))**2) )
                        if (distance1f <= 5 .or. distance2f <= 5) then
!                                write(*,*) "ii= " , ii, "and distance1f is" , distance1f, "distance2f= " , distance2f
                                dotdipolXf=(((Xwf(ii+1)+Xwf(ii+2))/2.0_dp)-Xwf(ii))
                                dotdipolYf=(((Ywf(ii+1)+Ywf(ii+2))/2.0_dp)-Ywf(ii))
                                dotdipolZf=(((Zwf(ii+1)+Zwf(ii+2))/2.0_dp)-Zwf(ii))
                                lengthf=SQRT(dotdipolXf**2+dotdipolYf**2+dotdipolZf**2)
                                Xf=(dotdipolXi/lengthi)*(dotdipolXf/lengthf)
                                Yf=(dotdipolYi/lengthi)*(dotdipolYf/lengthf)
                                Zf=(dotdipolZi/lengthi)*(dotdipolZf/lengthf)
                                unitdipol=Xf+Yf+Zf
                                dipolf(jj)=unitdipol
!                               write(*,*) "XH2i= " ,Xwi(ii+2), "XH1i= " ,Xwi(ii+1)
!                               write(*,*) "YH2i= " ,Ywi(ii+2), "YH1i= " ,Ywi(ii+1)
!                               write(*,*) "ZH2i= " ,Zwi(ii+2), "ZH1i= " ,Zwi(ii+1)
!                               write(*,*) "XOwi= " ,Xwi(ii), "YOwi= " ,Ywi(ii), "ZOwi= ", Zwi(ii)                              
!                               write(*,*)  "ii= ",ii, "X= " , dotdipolX,  "Y= " , dotdipolY, "Z= " , dotdipolZ
!                              write(*,*)  "ii= ",ii, "dipolf " , dipolf(jj)
                                jj=jj+1
                        end if
                end if
        end if        
end do
!timefunction=10
jj=jj-1; kk=kk-1; ! we started froom jj and kk=1; we need to remove one to have the correct # of waters

!timefunction=SUM(dipolf)                        !FLAG
!write(*,*) "SUM(dipolf)" ,timefunction          !FLAG
!timefunction=SUM(dipoli)                        !FLAG
!write(*,*) "SUM(dipoli)" ,timefunction          !FLAG
if (jj == 0 .or. kk == 0) then
        timefunction=0.0_dp
else
        timefunction=(SUM(dipolf)/jj)/(SUM(dipoli)/kk)
end if
!write(*,*) " # of H2o at zero is  ", kk, "# of H2O at t is ", jj, " time function is " , timefunction !FLAG
end subroutine RotationalCorrelationFunction

!deallocate(waterID,Xw,Yw,Zw,AdsID,Xa,Ya,Za,stat = alloc_err)

!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINE FOR READING CONFIG FORMAT FILE !!!!!!!!!!!!!!!!!!!!!!
subroutine ReadConfig(filename,metalMM,metalQM,Ads,Water,waterID,Xwat,Ywat,Zwat,AdsID,Xads,Yads,Zads)
implicit none

integer  ,  parameter                           ::      dp = selected_real_kind(15, 307)

!local variable
integer                                         ::      ierror,  error_flag, alloc_err
integer                                         ::      ii
character(len=256)                              ::      skip ! The variable forskipping lines
character(len=256)                              ::      str_2, str_3, str_4 !The variables for skipping strings in a line
!integer,   parameter                            ::      metalMM=1229
!integer,   parameter                            ::      metalQM=51
!integer,   parameter                            ::      Ads=10
!integer,   parameter                            ::      Water=2200

!parameter types and definition
character(len=30)    ,                intent(in)                  ::      filename 
integer,                              intent(in)                  ::      metalMM,metalQM,Ads,Water
real(dp), dimension(Water),  intent(out)                          ::      Xwat,Ywat,Zwat
character(len=256),dimension(Water), intent(out)                  ::      waterID
real(dp), dimension(Ads),  intent(out)                            ::      Xads,Yads,Zads 
character(len=256),dimension(Ads), intent(out)                    ::      AdsID


!Read the energies in the file
!allocate(metalMMID(metalMM))
!allocate(XmMM(metalMM))
!allocate(YmMM(metalMM))
!allocate(ZmMM(metalMM))
!allocate(waterID(Water))
!allocate(Xwat(waterlines))
!allocate(Ywat(waterlines))
!allocate(Zwat(waterlines))

!allocate(metalQMID(metalQM))
!allocate(XmQM(metalQM))
!allocate(YmQM(metalQM))
!allocate(ZmQM(metalQM))

!allocate(AdsID(Ads))
!allocate(Xads(Ads))
!allocate(Yads(Ads))
!allocate(Zads(Ads))

error_flag = 0

open(unit=77777, file=filename, status='old', action='read', iostat = ierror)

! skip the header lines, 2 lines
do ii=1,2
   read(77777,'(A)', iostat = ierror) skip
        if (ierror /=0 ) then                  !!!!FALG
          write(*,*) " Problem with reading file: ", filename
          error_flag = 1
          call EXIT(0)
        end if
  ! write(*,*) skip
end do
! Read Metal MM atoms, those are 1229 Pt atoms
do ii=1,metalMM
        read(77777,*, iostat = ierror) skip
                if (ierror /=0 ) then                  !!!!FALG
                write(*,*) " Problem with reading file: ", filename
                error_flag = 1
                call EXIT(0)
                end if
        !write(*,*)  skip
end do

do ii=1,Water
        read(77777,*, iostat = ierror) waterID(ii), str_2, str_3, str_4
                if (ierror /=0 ) then                  !!!!FALG
                write(*,*) " Problem with reading file: ", filename
                error_flag = 1
                call EXIT(0)
                end if
!        write(*,*)  waterID(ii)
        read(str_2,*) Xwat(ii)
!        write(*,*) Xwat(ii)
        read(str_3,*) Ywat(ii)
!        write(*,*) Ywat(ii)
        read(str_4,*) Zwat(ii)
!        write(*,*) Zwat(ii)
end do

do ii=1,metalQM
        read(77777,*, iostat = ierror) skip
                if (ierror /=0 ) then                  !!!!FALG
                write(*,*) " Problem with reading file: ", filename
                error_flag = 1
                call EXIT(0)
                end if
        !write(*,*) skip
end do
do ii=1,Ads
        read(77777,*, iostat = ierror) AdsID(ii), str_2, str_3, str_4
                if (ierror /=0 ) then                  !!!!FALG
                write(*,*) " Problem with reading file: ", filename
                error_flag = 1
                call EXIT(0)
                end if
!        write(*,*)  AdsID(ii)
        read(str_2,*) Xads(ii)
!        write(*,*) Xads(ii)
        read(str_3,*) Yads(ii)
!        write(*,*) Yads(ii)
        read(str_4,*) Zads(ii)
!        write(*,*) Zads(ii)
end do


close(77777)

!deallocate(metalMMID, XmMM, YmMM, ZmMM,waterID, stat = alloc_err)
!deallocate(metalQMID, XmQM, YmQM, ZmQM,AdsID, stat = alloc_err)

end subroutine ReadConfig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! aEND OF SUBROUTINE energyValues !!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTNE FOR COUNTING lines IN A FILE!!!!!!!!!!!!!
subroutine countline(filename, nline)
implicit none

integer  ,  parameter                   ::      dp = selected_real_kind(15, 307)

!local variable
integer                                 ::      ierror,  error_flag

!parameter types and definition
character(len=30)    ,  intent(in)                  ::      filename
integer              ,  intent(out)                 ::      nline


error_flag = 0

!read the # of lines in the file 
nline=0
open(1000, file=filename, status='old', action='read', iostat = ierror)
do
 read (1000, *, end=10)
        if (ierror /=0 ) then                  !!!!FALG
          write(*,*) " Problem with counting lines in file:  " , filename
           error_flag = 1
          call EXIT(0)
        end if
 nline=nline+1
end do
10 close(1000)

end subroutine countline
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF SUBROUTINE COUNTLINE !!!!!!!!!!!!!!!!!!!!!


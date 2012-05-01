PROGRAM ising
!Testing remote push

INTEGER,PARAMETER:: iterLen = 0!1000 !Length time averaging magnetization over
INTEGER,PARAMETER:: equilLen = 1!3000!3000 !Length time before averaging
INTEGER,PARAMETER:: SideN = 50 !Length of the sides of the lattice
INTEGER,PARAMETER:: TotalN =SideN**2!Total sites in the lattice
CHARACTER(LEN=*),PARAMETER:: fileName = "averaging_b_1.dat"

INTEGER,PARAMETER:: numTemps = 1!400!000 !Number of temps to iterate over
REAL*8, PARAMETER:: startTemp=1.5
REAL*8, PARAMETER:: endTemp=1.5!1.78
REAL*8 :: tempStep = (endTemp-startTemp)/numTemps
REAL*8 :: temp  !Temperature to run simulation at


INTEGER,DIMENSION(TotalN):: LatSpin !Holds the spin of each site on lattice
INTEGER,DIMENSION(TotalN):: flipSpins !Holds -1 for spins to flip, 1 to not

REAL*8 :: curMag    !Current magnetization of system
REAL*8,DIMENSION(iterLen) :: magDat !All curMag over iterLen period
REAL*8 :: curEn     !Current energy of system
REAL*8,DIMENSION(iterLen) :: enDat

INTEGER :: flip !1 if cluster is flipped up, -1 if cluster is flipped to down
INTEGER :: initialSite !Location of site that cluster is built around

INTEGER:: tempIter !Iterator for temperature loop
INTEGER:: i,l !Iterator for main loop
REAL*8 :: randInit !Random number for picking site

!--------------------------------------------------------------------------------

CALL Random()
temp=startTemp

DO tempIter =1, numTemps

    !Initialize----------------------------------------------
    CALL initialize()

    !Run Equilibration------------------------------------------------
    DO i = 1,equilLen
        CALL RANDOM_NUMBER(randInit)
        initialSite =MODULO(CEILING(1000*(randInit)*TotalN),TotalN)+1
    
        flip = -LatSpin(initialSite)
        flipSpins = wolffCluster(initialSite)     !Build the cluster
        LatSpin = LatSpin*flipSpins             !Flip spins of the cluster

    END DO

    !Run Main---------------------------------------------------------
    curMag = REAL(SUM(LatSpin))/TotalN
    DO i = 1,iterLen
        CALL RANDOM_NUMBER(randInit)
        initialSite =MODULO(CEILING(1000*(randInit)*TotalN),TotalN)+1

    
        flip = -LatSpin(initialSite)
        flipSpins = wolffCluster(initialSite)     !Build the cluster
        LatSpin = LatSpin*flipSpins             !Flip spins of the cluster

        curMag = getCurMag(flipSpins,flip,curMag)
        curEn = getCurEn()

        CALL UpdateAverages(ABS(curMag),magDat,i)
        CALL UpdateAverages(curEn,enDat,i)
    END DO

    !Call WriteAverage(magDat,"mag.dat")
    !CALL WriteVariance(magDat,"sus.dat")
    !CALL WriteAverage(EnDat, "en.dat")
    !CALL WriteVariance(EnDat, "Cv.dat")
    !CALL writeAll(magDat,EnDat,fileName)
END DO


!Run Outputs-------------------------------------------------------


!Testing-------------------------------------------------

DO i = 0, SideN-1
    DO l=1, SideN
        WRITE(*,"(1x,i4)",advance="no") LatSpin(i*SideN+l)
    END DO
    WRITE(*,*)
END DO

!OPEN(unit=27,status='replace',file='mag.dat')
!    DO i = 1, iterLen
!        WRITE(27,*) ABS(magDat(i))
!    END DO
!CLOSE(27)

!Subroutines---------------------------------------------
CONTAINS

SUBROUTINE initialize()
    IMPLICIT NONE
    
    CALL assignInitialSpins(LatSpin)
    temp = temp+tempStep  !Comment me for one temp
    magDat = 0

END SUBROUTINE

!Sets up initial spin configuration, currently all up
SUBROUTINE assignInitialSpins(lattice)
    INTEGER, DIMENSION(:) :: lattice
    REAL*8 :: randX

    lattice =1
    DO i = 1,TotalN
        CALL RANDOM_NUMBER(randX)
        lattice(i) = lattice(i)-2*FLOOR(randX+.5)
    END DO

END SUBROUTINE

!Updates averages based on sites that were flipped
SUBROUTINE UpdateAverages(curVal, dat,it)
    INTEGER :: it
    REAL*8 :: curVal
    REAL*8,DIMENSION(:) :: dat

    dat(it) = curVal
END SUBROUTINE

!Gets all 4 neighbours of a given site
FUNCTION getNeighbours(site)
    INTEGER :: site
    INTEGER,DIMENSION(4) :: getNeighbours
    INTEGER,DIMENSION(4) :: neighbours

    neighbours(1) = MODULO(site + SideN,TotalN +1) + FLOOR(REAL(site+SideN-1)/TotalN)
    neighbours(2) = MODULO(site - SideN + TotalN,TotalN +1)+ CEILING(REAL(site-SideN)/TotalN)
    neighbours(3) = (site + 1) - SideN*(CEILING(REAL(site+1)/SideN)-CEILING(REAL(site)/SideN))
    neighbours(4) = (site - 1) - SideN*(CEILING(REAL(site-1)/SideN)-CEILING(REAL(site)/SideN))

    getNeighbours = neighbours

    return
END FUNCTION

FUNCTION getCurEn()
    IMPLICIT NONE
    INTEGER :: site
    REAL*8 :: getCurEn

    getCurEn=0

    DO site = 1,TotalN
        getCurEn = getCurEn-latSpin(site)*SUM(latSpin(getNeighbours(site)))
    END DO

    RETURN
END FUNCTION

FUNCTION getCurMag(flips,flip,prevMag)
    INTEGER,DIMENSION(:) :: flips   !Return value from buildCluster (-1 if flip)
    INTEGER :: numFlips !Total number of sites flipped
    INTEGER :: flip !1 if cluster flipped to up, -1 if flipped to down
    REAL*8 :: prevMag !Previous magnification
    REAL*8 :: getCurMag
    

    flips = -FLOOR(flips/2.0)
    numFlips = SUM(flips)

    getCurMag = prevMag + 2.0*flip*numFlips/TotalN
    return 

END FUNCTION

SUBROUTINE random

    IMPLICIT NONE

    INTEGER :: seed, i
    REAL*8 :: x

    OPEN(UNIT=26,file='seed')
        READ(26,*) seed
    CLOSE(26)

    DO i=1,seed
        CALL random_number(x)
    END DO

    OPEN(unit=26,status='replace',file='seed')
        WRITE(26,*) CEILING(1000*x)
    CLOSE(26)

END SUBROUTINE

SUBROUTINE writeAll(mag,en,fName)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: fName
    REAL*8, DIMENSION(:) :: mag
    REAL*8, DIMENSION(:) :: en
    REAL*8 :: avgMag, avgMag2, varMag
    REAL*8 :: avgEn, avgEn2, varEn

    avgMag2 = SUM(mag**2)/iterLen
    avgMag = SUM(mag)/iterLen
    varMag = avgMag2-avgMag**2

    avgEn2 = SUM(en**2)/iterLen
    avgEn = SUM(en)/iterLen
    varEn = avgEn2-avgEn**2

    OPEN(UNIT=26,file=fName,POSITION="APPEND")
        WRITE(26,*) temp,avgMag, varMag, avgEn, varEn
    CLOSE(26)

END SUBROUTINE

SUBROUTINE WriteVariance(dataArray,fName)
    IMPLICIT NONE
    REAL *8, DIMENSION(:) :: dataArray
    CHARACTER(LEN=*) :: fName
    REAL *8 :: avg, avg2
    REAL *8 :: variance

    avg2 = SUM(dataArray**2)/iterLen
    avg = SUM(dataArray)/iterLen
    variance = avg2-avg**2


    OPEN(UNIT=26,file=fName,POSITION="APPEND")
        WRITE(26,*) temp, variance
    CLOSE(26)

END SUBROUTINE

SUBROUTINE WriteAverage(dataArray,fName)
    IMPLICIT NONE
    REAL*8,DIMENSION(:) :: dataArray
    CHARACTER(LEN=*) :: fName

    OPEN(UNIT=26,file=fName,POSITION="APPEND")
        WRITE(26,*) temp,  SUM(dataArray)/iterLen
    CLOSE(26)
END SUBROUTINE

SUBROUTINE TemperatureInput(temperature)
    IMPLICIT NONE
    REAL*8 :: temperature
    
    READ(*,*) temperature

END SUBROUTINE

!WolffCluster Subroutines and functions-----------------------------------------

!Builds a cluster using the wolff algorithm
FUNCTION wolffCluster( initSite)
    INTEGER :: j    !Integer to iterate through queue (head)

    INTEGER :: initSpin !Spin sites will flip from (spin of init before flip)
    INTEGER :: initSite !Initial site to flip
    INTEGER :: curSite  !current site being added

    INTEGER,DIMENSION(totalN) :: inQueue    !0 if site not in queue, 1 if it is
    INTEGER,DIMENSION(totalN) :: wolffCluster   !Current cluster of spins 

    INTEGER,DIMENSION(totalN*5) :: Queue    !Queue of sites to try adding
    INTEGER :: QueuePos     !Site of last added queue element (tail)

    !Initializing Queue
    Queue = 0
    Queue(1:4) = getNeighbours(initSite)
    QueuePos = 5

    !Initializing status arrays
    inQueue = 0
    stat=0
    wolffCluster = 1

    initSpin = latSpin(initSite)
    wolffCluster(initSite) = -1

    !Main loop to build cluster
    j=1


    DO WHILE (j.LT.QueuePos)
        curSite = Queue(j)
        IF (wolffCluster(curSite).EQ.(1)) THEN
            CALL tryAddSite(curSite, Queue,QueuePos,wolffCluster,inQueue, initSpin)
        END IF
        j=j+1
    END DO


    return  !Returns cluster -1 to flip, 1 to not flip
END FUNCTION


!Tries to add site to cluster
SUBROUTINE tryAddSite(site, qSites,qPos,wCluster,inQ,iSpin)
    IMPLICIT NONE
    INTEGER :: iSpin !Spin sites will flip from (spin of initSite before flip)
    INTEGER :: site !Current site to try and add
    INTEGER :: qPos !Position of last element added to queue (tail)
    INTEGER :: k    !Loop var to iterate through neighbours
    INTEGER, DIMENSION(4) :: neighbours !Array of neighbours of site
    INTEGER, DIMENSION(:) :: qsites !Currentqueue
    INTEGER, DIMENSION(:) :: inQ    !Array 0 if site not in queue yet
    INTEGER, DIMENSION(:) :: wCluster   !Current cluster
    REAL*8 :: randX
    REAL*8 :: prob
    REAL*8 :: J = 1.0d0

    CALL RANDOM_NUMBER(randX)

    neighbours = getNeighbours(site)    
    prob = 1.0-EXP(-2.0d0*J/temp)

    IF ((LatSpin(site).EQ.iSpin).AND.(randX.LT.prob)) THEN   

        wCluster(site) = -2     !Adds site to cluster

        
        !Loop to add neighbours to queue
        DO k=1,4
            IF (inQ(neighbours(k)).LT.4) THEN   !Only add if not in queue
                qSites(qPos) = Neighbours(k)
                qPos = qPos +1
                inQ(Neighbours(k)) = inQ(Neighbours(k))+1
            END IF
        END DO

    END IF

END SUBROUTINE

END PROGRAM



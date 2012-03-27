PROGRAM ising

INTEGER,PARAMETER:: SideN=10 !Length of the sides of the lattice
INTEGER,PARAMETER:: TotalN = SideN**2 !Total sites in the lattice

INTEGER,DIMENSION(TotalN):: LatSpin !Holds the spin of each site on lattice
INTEGER,DIMENSION(TotalN):: flipSpins !Holds -1 for spins to flip, 1 to not

REAL*8 :: curMag    !Current magnetization of system

INTEGER :: flip !1 if cluster is flipped up, -1 if cluster is flipped to down
INTEGER :: initialSite !Location of site that cluster is built around

INTEGER:: i,l !Iterator for main loop
REAL :: temp=5  !Temperature to run simulation at


!Initialize----------------------------------------------
CALL assignInitialSpins(LatSpin)

curMag = 1.0

!Run Main------------------------------------------------
DO i = 1,1
    initialSite = 1

    
    flip = LatSpin(initialSite)
    flipSpins = wolffCluster(initialSite)     !Build the cluster
    LatSpin = LatSpin*flipSpins             !Flip spins of the cluster
  


    curMag = getCurMag(flipSpins,flip,curMag)      
    CALL UpdateAverages(curMag)          !Update the averages based on flips


END DO


!Testing-------------------------------------------------
DO i = 0, SideN-1
    DO l=1, SideN
        WRITE(*,"(1x,i4)",advance="no") LatSpin(i*4+l)
    END DO
    WRITE(*,*)
END DO

!WRITE(*,*) curMag


!Subroutines---------------------------------------------
CONTAINS

!Sets up initial spin configuration, currently all up
SUBROUTINE assignInitialSpins(lattice)
    INTEGER, DIMENSION(:) :: lattice
    lattice = 1

END SUBROUTINE

!Builds a cluster using the wolff algorithm
FUNCTION wolffCluster( initSite)
    INTEGER :: j    !Integer to iterate through queue (head)

    INTEGER :: initSpin !Spin sites will flip from (spin of init before flip)
    INTEGER :: initSite !Initial site to flip
    INTEGER :: curSite  !current site being added

    INTEGER,DIMENSION(totalN) :: inQueue    !0 if site not in queue, 1 if it is
    INTEGER,DIMENSION(totalN) :: wolffCluster   !Current cluster of spins 

    INTEGER,DIMENSION(totalN*4) :: Queue    !Queue of sites to try adding
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

    !DO j=1,1
    !DO WHILE (j.LT.QueuePos-1)
         
        !j=j+1

        curSite = Queue(j)
        CALL tryAddSite(curSite, Queue,QueuePos,wolffCluster,inQueue, initSpin)

    !END DO

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
    INTEGER, DIMENSION(:) :: qsites !Current queue
    INTEGER, DIMENSION(:) :: inQ    !Array 0 if site not in queue yet
    INTEGER, DIMENSION(:) :: wCluster   !Current cluster
    REAL*8 :: randX
    REAL*8 :: prob
    REAL*8 :: J = 1.0d0

    CALL RANDOM_NUMBER(randX)
    prob =1 !1.0-EXP(-2.0d0*J/temp)

    IF ((LatSpin(site).EQ.iSpin).AND.(randX.LT.prob)) THEN   

        wCluster(site) = -1     !Adds site to cluster
        neighbours = getNeighbours(site)    
        
        !Loop to add neighbours to queue
        DO k=1,4
            IF (inQ(neighbours(k)).EQ.0) THEN   !Only add if not in queue
                qSites(qPos) = Neighbours(k)
                qPos = qPos +1
                inQ(Neighbours(k)) = 1
            END IF
        END DO

    END IF

END SUBROUTINE

!Updates averages based on sites that were flipped
SUBROUTINE UpdateAverages(curM)
    REAL*8 :: curM

END SUBROUTINE

!Gets all 4 neighbours of a given site
FUNCTION getNeighbours(site)
    INTEGER :: site
    INTEGER,DIMENSION(4) :: getNeighbours
    INTEGER,DIMENSION(4) :: neighbours

    neighbours(1) = MODULO(site + SideN,TotalN +1) + FLOOR(REAL(site+SideN-1)/TotalN)
    neighbours(2) = MODULO(site - SideN + TotalN,TotalN +1)+ CEILING(REAL(site-SideN)/TotalN)
    neighbours(3) = (site + 1) - 4*(CEILING(REAL(site+1)/SideN)-CEILING(REAL(site)/SideN))
    neighbours(4) = (site - 1) - 4*(CEILING(REAL(site-1)/SideN)-CEILING(REAL(site)/SideN))

    getNeighbours = neighbours

    return
END FUNCTION

FUNCTION getCurMag(flips,flip,prevMag)
    INTEGER,DIMENSION(:) :: flips   !Return value from buildCluster (-1 if flip)
    INTEGER :: numFlips !Total number of sites flipped
    INTEGER :: flip !1 if cluster flipped to up, -1 if flipped to down
    REAL*8 :: prevMag !Previous magnification
    REAL*8 :: getCurMag
    
    
    flips = FLOOR(flips/2.0)
    numFlips = SUM(flips)

    getCurMag = prevMag + 2.0*flip*numFlips/REAL(TotalN)
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

END PROGRAM



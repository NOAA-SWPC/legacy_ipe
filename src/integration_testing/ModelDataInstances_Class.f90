! ModelDataInstances_Class.f90
! 
!  Copyright (2017), Joseph Schoonover, Cooperative Institute for Research in
!  Environmental Sciences, NOAA, (joseph.schoonover@noaa.gov)
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE ModelDataInstances_Class

!src/common/
USE Module_Precision

IMPLICIT NONE

   INTEGER, PARAMETER      :: strLen = 30
   CHARACTER(6), PARAMETER :: strFMT = "(A30)"

   TYPE ModelDataInstance
      CHARACTER(strLen)               :: moduleName
      CHARACTER(strLen)               :: subroutineName
      CHARACTER(strLen)               :: statusCheckName
      INTEGER                         :: lineNumber
      INTEGER                         :: instanceID
      INTEGER                         :: nObs
      INTEGER                         :: arraySize
      REAL(real_prec), ALLOCATABLE    :: array(:) 
      TYPE(ModelDataInstance),POINTER :: next

   END TYPE ModelDataInstance

   TYPE ModelDataInstances
      TYPE(ModelDataInstance), POINTER :: head, current, tail

      CONTAINS

      PROCEDURE :: Build       => Build_ModelDataInstances
      PROCEDURE :: Trash       => Trash_ModelDataInstances 
      PROCEDURE :: AddInstance => AddInstance_ModelDataInstances 

      PROCEDURE :: SetNames            => SetNames_ModelDataInstances
      PROCEDURE :: GetNames            => GetNames_ModelDataInstances
      PROCEDURE :: PointToInstance     => PointToInstance_ModelDataInstances
      PROCEDURE :: ThereAreNoInstances
        
      PROCEDURE :: Write_ModelDataInstances
      PROCEDURE :: Read_ModelDataInstances

    END TYPE ModelDataInstances

 INTEGER, PRIVATE, PARAMETER :: defaultNameLength = 40
 CONTAINS
!
!
!==================================================================================================!
!---------------------------- CONSTRUCTOR/DESTRUCTOR ROUTINES -------------------------------------!
!==================================================================================================!
!
 SUBROUTINE Build_ModelDataInstances( theInstances  )

   IMPLICIT NONE
   CLASS( ModelDataInstances ) :: theInstances
   
     theInstances % head => NULL( )
     ! Point the tail to null
     theInstances % tail => NULL()
     ! Set the current position to Null
     theInstances % current => NULL( )


 END SUBROUTINE Build_ModelDataInstances
!
 SUBROUTINE Trash_ModelDataInstances( theInstances )

   IMPLICIT NONE
   CLASS( ModelDataInstances ) :: theInstances
   ! LOCAL 
   TYPE( ModelDataInstance ), POINTER :: pNext
   
     ! Set the current position of the list to the head
     theInstances % current => theInstances % head
     
     ! Scroll through the list until the current position is nullified
     DO WHILE ( ASSOCIATED( theInstances % current ) )

        ! temporarily point to the next in the list
        pNext => theInstances % current % next 

        ! Deallocate memory pointed to by current position
        DEALLOCATE( theInstances % current % array )
        DEALLOCATE( theInstances % current ) 

        ! Update current position
        theInstances % current => pNext 

     ENDDO
      

 END SUBROUTINE Trash_ModelDataInstances
!
!
!==================================================================================================!
!---------------------------------- ACCESSOR ROUTINES ---------------------------------------------!
!==================================================================================================!
!
! 
 SUBROUTINE SetNames_ModelDataInstances( theInstances, moduleName, subroutineName, statusCheckName, instanceID )

   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(inout) :: theInstances
   CHARACTER(*), INTENT(IN)                   :: moduleName, subroutineName, statusCheckName
   INTEGER, OPTIONAL, INTENT(IN)              :: instanceID
      
      IF( PRESENT( instanceID ) )THEN
         CALL theInstances % PointToInstance( instanceID )  
      ENDIF
      
      theInstances % current % moduleName      = TRIM(moduleName)
      theInstances % current % subroutineName  = TRIM(subroutineName)
      theInstances % current % statusCheckName = TRIM(statusCheckName)

 END SUBROUTINE SetNames_ModelDataInstances
!
 SUBROUTINE GetNames_ModelDataInstances( theInstances, moduleName, subroutineName, statusCheckName, instanceID )

   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(in) :: theInstances
   CHARACTER(*), INTENT(out)               :: moduleName, subroutineName, statusCheckName
   INTEGER, OPTIONAL, INTENT(IN)           :: instanceID
      
      IF( PRESENT( instanceID ) )THEN
         CALL theInstances % PointToInstance( instanceID )
      ENDIF
      
      moduleName      = theInstances % current % moduleName
      subroutineName  = theInstances % current % subroutineName
      statusCheckName = theInstances % current % statusCheckName

 END SUBROUTINE GetNames_ModelDataInstances
!
!
!==================================================================================================!
!-------------------------------- Linked-List Type Operations -------------------------------------!
!==================================================================================================!
!
! 
 FUNCTION ThereAreNoInstances( theInstances ) RESULT( TorF )
  IMPLICIT NONE
  CLASS( ModelDataInstances ) :: theInstances
  LOGICAL              :: TorF

     TorF = .NOT.( ASSOCIATED( theInstances % head  ) )
     
 END FUNCTION ThereAreNoInstances
!
 SUBROUTINE AddInstance_ModelDataInstances( theInstances, &
                                            moduleName, &
                                            subroutineName, &
                                            statusCheckName, &
                                            lineNumber,& 
                                            arraySize )
 
   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(inout) :: theInstances
   CHARACTER(*)                               :: moduleName
   CHARACTER(*)                               :: subroutineName
   CHARACTER(*)                               :: statusCheckName
   INTEGER                                    :: lineNumber, arraySize
   ! LOCAL
   INTEGER :: allocationStatus

     ! Check to see if this list is empty
     IF( theInstances % ThereAreNoInstances() )THEN
     
        ALLOCATE( theInstances % head, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE ModelDataInstances_Class.f90 : S/R AddInstance : Memory not allocated for next entry in list.'
        ENDIF      
      
        ! Point the current position to the head
        theInstances % current => theInstances % head
        ! Set the data
        CALL theInstances % SetNames( moduleName, subroutineName, statusCheckName  )
        
        theInstances % current % instanceID = CharToIntHashFunction( theInstances % current % statusCheckName )
        theInstances % current % nObs       = 0
        theInstances % current % arraySize   = arraySize
        ALLOCATE( theInstances % current % array(1:arraySize) )
        theInstances % current % array = 0.0_real_prec
        
        ! Point the next to null and the tail to current
        theInstances % current % next => NULL( )
        theInstances % tail => theInstances % current
        
     ELSE ! the list is not empty
    
        ! Then we allocate space for the next item in the list    
        ALLOCATE( theInstances % tail % next, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE ModelDataInstances_Class.f90 : S/R AddInstance : Memory not allocated for next entry in list.'
        ENDIF      
        
        ! Reassign the tail
        theInstances % tail => theInstances % tail % next
        
        ! Set the current to the tail
        theInstances % current => theInstances % tail
  
        ! Fill in the data
        CALL theInstances % SetNames( moduleName, subroutineName, statusCheckName  )
        
        ! Fill in the key information
        theInstances % current % instanceID = CharToIntHashFunction( theInstances % current % statusCheckName )
        theInstances % current % nObs       = 0
        theInstances % current % arraySize   = arraySize
        ALLOCATE( theInstances % current % array(1:arraySize) )
        theInstances % current % array = 0.0_real_prec
        
        ! Point the next to null and the tail to current
        theInstances % current % next => NULL( )
        
     ENDIF

 END SUBROUTINE AddInstance_ModelDataInstances
!
 SUBROUTINE PointToInstance_ModelDataInstances( theInstances, instanceID, instanceFound )
 
   IMPLICIT NONE
   CLASS( ModelDataInstances )    :: theInstances
   INTEGER, INTENT(IN)            :: instanceID
   LOGICAL, OPTIONAL, INTENT(OUT) :: instanceFound
   
   
      theInstances % current => theInstances % head ! Point to the head of the list
      IF( PRESENT( instanceFound ) )THEN
        instanceFound = .FALSE.
      ENDIF
      
      DO WHILE(ASSOCIATED(theInstances % current))
      
         IF( theInstances % current % instanceID == instanceID )THEN
            IF( PRESENT( instanceFound ) )THEN
              instanceFound = .TRUE.
            ENDIF
            EXIT
         ENDIF
         theInstances % current => theInstances % current % next
       
      ENDDO
      
      
 END SUBROUTINE PointToInstance_ModelDataInstances
!
 SUBROUTINE Write_ModelDataInstances( theInstances, baseFileName )

   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(INOUT) :: theInstances 
   CHARACTER(200), INTENT(IN)                 :: baseFileName
   ! LOCAL
   INTEGER :: k, fUnit, fUnit2, recID, i
   CHARACTER(3) :: countChar 

      theInstances % current => theInstances % head
      WRITE( countChar, '(I3.3)' ) theInstances % current % nObs
      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(baseFileName)//'.instance.header', &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            ACTION = 'WRITE' )
            
      OPEN( UNIT = NewUnit(fUnit2), &
            FILE = TRIM(baseFileName)//'.'//countChar//'.instance', &
            FORM = 'UNFORMATTED', &
            ACCESS = 'DIRECT', &
            ACTION = 'WRITE', &
            RECL   = real_prec ) 

      k     = 0
      recID = 0
      DO WHILE( ASSOCIATED(theInstances % current) )
         k = k+1


         WRITE(fUnit,*) TRIM( theInstances % current % moduleName )
         WRITE(fUnit,*) TRIM( theInstances % current % subroutineName )
         WRITE(fUnit,*) TRIM( theInstances % current % statusCheckName )
         WRITE(fUnit,*) theInstances % current % lineNumber
         WRITE(fUnit,*) theInstances % current % arraySize
         WRITE(fUnit,*) theInstances % current % instanceID
         WRITE(fUnit,*) '------------------------------------------------------------'

         theInstances % current => theInstances % current % next

       
         DO i = 1, theInstances % current % arraySize
            recID = recID + 1
            WRITE( fUnit2, REC=recID ) theInstances % current  % array(i)
         ENDDO

      ENDDO
      CLOSE(fUnit)
      CLOSE(fUnit2)


 END SUBROUTINE Write_ModelDataInstances
!
 SUBROUTINE Read_ModelDataInstances( theInstances, baseFileName )

   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(INOUT) :: theInstances 
   CHARACTER(200), INTENT(IN)                 :: baseFileName
   ! LOCAL
   INTEGER        :: k, fUnit, fUnit2, recID, i
   CHARACTER(3)   :: countChar 
   CHARACTER(strLen) :: moduleName, subroutineName, statusCheckName, dummyChar
   INTEGER        :: lineNumber, arraySize, instanceID

      WRITE( countChar, '(I3.3)' ) theInstances % current % nObs
      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(baseFileName)//'.instance.header', &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            ACTION = 'READ' )
            
      OPEN( UNIT = NewUnit(fUnit2), &
            FILE = TRIM(baseFileName)//'.'//countChar//'.instance', &
            FORM = 'UNFORMATTED', &
            ACCESS = 'DIRECT', &
            ACTION = 'READ', &
            RECL   = real_prec ) 

      k     = 0
      recID = 0
      DO WHILE( ASSOCIATED(theInstances % current) )
         k = k+1


         READ(fUnit,strFMT) moduleName 
         READ(fUnit,strFMT) subroutineName
         READ(fUnit,strFMT) statusCheckName
         READ(fUnit,*) lineNumber
         READ(fUnit,*) arraySize
         READ(fUnit,*) instanceID
         READ(fUnit,strFMT) dummyChar

         CALL theInstances % AddInstance( moduleName, &
                                          subroutineName, &
                                          statusCheckName, &
                                          lineNumber, &
                                          arraySize )
                                          

       
         DO i = 1, theInstances % current % arraySize
            recID = recID + 1
            READ( fUnit2, REC=recID ) theInstances % current  % array(i)
         ENDDO

      ENDDO
      CLOSE(fUnit)
      CLOSE(fUnit2)


 END SUBROUTINE Read_ModelDataInstances
!
 FUNCTION CharToIntHashFunction( inputChar ) RESULT( hash )
    IMPLICIT NONE
    CHARACTER(strLen) :: inputChar
    INTEGER           :: hash
    ! Local
    CHARACTER(strLen) :: localChar
    INTEGER           :: i

       localChar = UpperCase( inputChar )
      
       hash = 5381
       DO i = 1, LEN_TRIM(localChar)
          hash = ( ISHFT(hash,5)+hash ) + ICHAR( localChar(i:i) )
       ENDDO

 END FUNCTION CharToIntHashFunction
!
 FUNCTION UpperCase( str ) RESULT( upper )

    IMPLICIT NONE
    CHARACTER(strLen), INTENT(IN) :: str
    CHARACTER(strLen)           :: Upper

    INTEGER :: ic, i

    CHARACTER(27), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '
    CHARACTER(27), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz '

    DO i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        IF (ic > 0)THEN
         Upper(i:i) = cap(ic:ic)
        ELSE
         Upper(i:i) = str(i:i)
        ENDIF
    ENDDO

 END FUNCTION UpperCase
!
 INTEGER FUNCTION NewUnit(thisunit)
 
   IMPLICIT NONE
   INTEGER, INTENT(out), OPTIONAL :: thisunit
   ! Local
   INTEGER, PARAMETER :: unitMin=100, unitMax=1000
   LOGICAL :: isopened
   INTEGER :: iUnit

     newunit=-1

     DO iUnit=unitMin, unitMax
        ! Check to see IF this UNIT is opened
        INQUIRE(UNIT=iUnit,opened=isopened)
        IF( .not. isopened )THEN
           newunit=iUnit
           EXIT
        ENDIF
     ENDDO

     IF (PRESENT(thisunit)) thisunit = newunit
 
END FUNCTION NewUnit
END MODULE ModelDataInstances_Class

! ##############################################################################################

! @@@@   MODULES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################

module Global
implicit none

integer :: nInd, nSnp, nIndLH, maxSibSize, MaxMismatch, maxOppHom, &
nRounds, nC(2), nAgeClasses, nPairs
integer, allocatable, dimension(:) :: Sex, BY, PairType, nFS
integer, allocatable, dimension(:,:) :: Genos, AgeDiff, Parent, OppHomM, &
    nS, PairID, FSID
integer, allocatable, dimension(:,:,:) :: SibID, GpID
double precision :: thLR, thLRrel, Er
double precision, allocatable, dimension(:) ::  Lind, PairDLLR, AF
double precision, allocatable, dimension(:,:) :: AHWE, OHWE, LLR_O, LindX, &
    LR_parent, AgePriorM, CLL
double precision, allocatable, dimension(:,:,:) :: AKAP, OKOP,AKOP,OKAP, &
AcO, LR_GP, LindG, PHS, PFS
double precision, allocatable, dimension(:,:,:,:) :: OKO2P, OKOAP,AKO2P,AKOAP,OKA2P, DumP
double precision, allocatable, dimension(:,:,:,:,:) :: XPr
double precision :: AKA2P(3,3,3)
 character(len=30), allocatable, dimension(:) :: Id, NameLH
 character(len=200) :: ParentageFileName, PedigreeFileName, AgePriorFileName
 character(len=2) :: DumPrefix(2)
 character(len=30), allocatable, dimension(:,:) :: ParentName
 
  contains
pure function MaxLL(V)
double precision, intent(IN) :: V(:)
double precision :: MaxLL

if (ANY(V < 0)) then
    MaxLL = MAXVAL(V, mask = V<0, DIM=1)
else
    MaxLL = MINVAL(V, DIM=1)  ! 777: can't do; 888: already is; 999: not calc'd
endif
end function MaxLL

end module Global

! ##############################################################################################

! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
! Made F conformant by Walt Brainerd

! Adapted by J Huisman (j.huisman@ed.ac.uk) to output rank, in order to enable 
! sorting of parallel vectors, and changed to decreasing rather than increasing order

module qsort_c_module
implicit none
public :: QsortC
private :: Partition

 contains
recursive subroutine QsortC(A, Rank)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: Rank
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq, Rank)
     call QsortC(A(:iq-1), Rank(:iq-1))
     call QsortC(A(iq:), Rank(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker, Rank)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: Rank
  integer, intent(out) :: marker
  integer :: i, j, TmpI
  double precision :: temp
  double precision :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp 
        
        TmpI = Rank(i) 
        Rank(i) = Rank(j)
        Rank(j) = TmpI
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module


! ##############################################################################################

! @@@@   PROGRAMS   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################


subroutine duplicates(nDupGenoID, nDupLhID, nDupGenos, nSexless)
!program Duplicates
use Global
implicit none

integer, intent(OUT) :: nDupGenoID, nDupLhID, nDupGenos, nSexless   
integer :: i, j, l, Match, CountMismatch
integer, allocatable, dimension(:,:) :: dupGenoIDs, dupLhIDs, DupGenos
integer, allocatable, dimension(:) :: Sexless, nMisMatch

 call ReadData
 
nDupGenoID = 0
nDupLhID = 0
nDupGenos = 0
nSexless = 0
allocate(DupGenoIDs(nInd,2))
allocate(DupLhIDs(nIndLH,2))
allocate(DupGenos(nInd,2))
allocate(Sexless(nInd))
allocate(nMismatch(nInd))

! subroutine intpr (label, nchar, data, ndata)
! call intpr ("nInd Geno: ", 11, nInd, 1)
! call intpr ("nInd LH: ", 9, nIndLH, 1)
! any duplicate IDs?
do i=1,nInd-1
    do j=i+1, nInd
        if (Id(i) == Id(j)) then
            nDupGenoID = nDupGenoID + 1
            dupGenoIDs(nDupGenoID,1) = i
            dupGenoIDs(nDupGenoID,2) = j
        endif
    enddo
enddo
 
do i=1,nIndLH-1
    do j=i+1, nIndLH
        if (NameLH(i) == NameLH(j)) then
            nDupLhID = nDupLhID + 1
            DupLhIDs(nDupLhID,1) = i
            DupLhIDs(nDupLhID,2) = j
        endif
    enddo
enddo

! identical genotypes?
do i=1,nInd-1
    do j=i+1, nInd
        Match=1
        CountMismatch=0
        do l=1, nSnp
            if (Genos(l,i)==-9 .or. Genos(l,j)==-9) cycle
            if (Genos(l,i) /= Genos(l,j)) then
                CountMismatch=CountMismatch+1
                if (CountMismatch > MaxMismatch) then
                    Match=0
                    exit
                endif
            endif
        enddo
        if (Match==1) then
            nDupGenos = nDupGenos + 1
            DupGenos(nDupGenos,1) = i
            DupGenos(nDupGenos,2) = j
            nMisMatch(nDupGenos) = CountMismatch
        endif
    enddo
enddo

! check if all genotyped individuals have a gender in the lifehistory data
do i=1,nInd
    if (Sex(i)/=1 .and. Sex(i)/=2) then
        nSexless = nSexless + 1
        Sexless(nSexless) = i
    endif
enddo

! call intpr ( "dup: ",5, (/ nDupGenoID, nDupLhID, nDupGenos, nSexless /), 4)
! ##########################

open (unit=201,file="DuplicatesFound.txt",status="replace")
write (201, '(a15, a10, a30, a10, a30, a10)') "Type", "Row1", "ID1", "Row2", "ID2", "nDiffer"
if (nDupGenoID>0) then
    do i=1,nDupGenoID
        write (201,'(a15, i10, " ", a30, i10, " ",a30)') "GenoID", dupGenoIDs(i,1), Id(dupGenoIDs(i,1)), &
        dupGenoIDs(i,2), Id(dupGenoIDs(i,2))
    enddo   
endif
if (nDupLhID>0) then
    do i=1,nDupLhID
        write (201,'(a15, i10," ", a30, i10, " ",a30)') "LifehistID", DupLhIDs(i,1), NameLH(DupLhIDs(i,1)), &
        DupLhIDs(i,2), NameLH(DupLhIDs(i,2))
    enddo   
endif
if (nDupGenos>0) then
    do i=1,nDupGenos
        write (201,'(a15, i10, " ",a30, i10, " ",a30, i10)') "Genotype", dupGenos(i,1), Id(dupGenos(i,1)), &
        dupGenos(i,2), Id(dupGenos(i,2)), nMismatch(i)  
    enddo   
endif
if (nSexless>0) then
    do i=1, nSexless
        write (201,'(a15, i10," ", a30, i10, " ",a30)') "NoSex", Sexless(i), Id(Sexless(i)), 0, "NA"
    enddo
endif
 close (201)

if (allocated (DupGenoIDs)) deallocate(DupGenoIDs)
if (allocated (DupLhIDs))  deallocate(DupLhIDs)
if (allocated (DupGenos))  deallocate(DupGenos)
if (allocated (Sexless))   deallocate(Sexless)
if (allocated (nMismatch)) deallocate(nMismatch)

 call DeAllocAll

end subroutine duplicates
!end program Duplicates

! ##############################################################################################

subroutine parents(nParents, nAmbiguous)
!program AssignParents
use qsort_c_module
use Global
implicit none

integer, intent(INOUT) :: nParents(2), nAmbiguous
integer :: i, j, k, Round, isP(2), topX, maybe
integer, allocatable, dimension(:) :: BYRank
integer, allocatable, dimension(:,:) :: OppHomDF
double precision, allocatable, dimension(:) :: SortBY
double precision :: LLtmp(7,3), dLL
 character(len=2) :: RelName(8)
               
 call ReadData
 call PrecalcProbs

maxOppHom = MaxMismatch - FLOOR(-nSNP * Er)   ! round up to nearest integer  

allocate(OppHomM(nInd, nInd))
OppHomM = 999 
allocate(LLR_O(nInd, nInd))
LLR_O = 999
allocate(Parent(nInd,2)) 
Parent = 0
allocate(OppHomDF(nInd,2))
OppHomDF = -9
allocate(SortBY(nInd))
allocate(BYRank(nInd))

nParents = 0
nAmbiguous = 0
nC = 0

!============================

do i=1,nInd
    call CalcLind(i)
enddo
 call intpr ( "Parentage ... ", -1, 0, 0)
 call dblepr("Initial total LL : ", -1, SUM(Lind), 1) 

 call CalcOppHom   ! also checks no. SNPs typed in both
 
do i=1, nInd-1
    do j=i+1,nInd 
        if (OppHomM(i,j) > maxOppHom) cycle    
        if (AgeDiff(i,j) /= 999) then
            if (AgeDiff(i,j) > 0) then  ! j older than i
                call CalcPO(i, j, LLR_O(i,j))  ! LLR PO/U
            else if (AgeDiff(i,j) < 0) then
                call CalcPO(j, i, LLR_O(j,i))
            endif
        else
            call CalcPO(i, j, LLR_O(i,j))
            call CalcPO(j, i, LLR_O(j,i))
        endif  
    enddo
enddo

! get birthyear ranking (increasing)
 SortBY = REAL(BY, 8)
 WHERE (SortBY < 0) SortBY = HUGE(0.0D0) 
 BYRank = (/ (i, i=1, nInd, 1) /)
 call QsortC(SortBy, BYRank)
 
! loop through, assign parents where possible. start with oldest.
allocate(LR_parent(nInd,3))
LR_parent = 999

do Round=1,3
    call Parentage(BYrank)
     
    do i=1,nInd
        call CalcLind(i)
    enddo
    
    ! assign sex if assigned as parent >= 2 times
    do i=1,nInd
        if (Sex(i)==3) then
            isP = 0
            do k=1,2
                do j=1,nInd
                    if (Parent(j,k) == i) then
                        isP(k) = isP(k) + 1
                    endif
                enddo
            enddo
            if (isP(1)>0 .and. isP(2)>0) then
                call rwarn("Assigned as both dam & sire: "//ID(i))
            else
                do k=1,2
                    if (isP(k)>1) then
                        Sex(i) = k
                    endif
                enddo
            endif
        endif
    enddo
enddo

 call UpdateAllProbs
 call dblepr("Post-parentage total LL : ", -1, SUM(Lind), 1) 
 call CalcParentLLR

do i=1,nInd
    do k=1,2
        if (Parent(i,k)>0) then
            ParentName(i,k) = Id(Parent(i,k))
            nParents(k) = nParents(k) + 1
            OppHomDF(i,k) = OppHomM(i, Parent(i,k))  
        endif
    enddo
enddo
    
open (unit=201,file=trim(ParentageFileName), status="unknown")  ! "Parents_assigned.txt"
write (201, '(3a30, 5a12, 3a6)') "ID", "Dam", "Sire", "LLR_dam", "LLR_sire", "LLR_pair", & 
"OH_Mother", "OH_Father", "RowO", "RowD", "RowS"   
do i=1,nInd
    write (201,'(3a30, 3f10.2, 2i12, 3i6)') Id(i), ParentName(i,1:2), LR_parent(i,1:3), &
    OppHomDF(i,1:2), i, Parent(i,1:2)
enddo   
 close (201)
! last 3 columns (9-11) read in for sibship assignment - update subroutine ReadParents when changing number of columns!!

! output non-assigned likely parents-offspring pairs. (unknown sex or birth year etc.) 
RelName = (/ "PO", "FS", "HS", "GP", "FA", "HA", "U ", "XX" /) 
open (unit=201,file="Unassigned_relatives_par.txt",status="unknown")
write (201, '(2a15, 2a5, 4a10)') "ID1", "ID2", "Sex1", "Sex2", "AgeDif", "BestRel", "LLR_R_U", "LLR_R1_R2"
do i=1,nInd-1
    do j=i+1,nInd
        if (OppHomM(i,j) > maxOppHom .or. OppHomM(i,j)<0) cycle
        if (LLR_O(i,j)==999 .or. LLR_O(i,j)< 0) cycle
        if (ANY(Parent(i,:)==j) .or. ANY(Parent(j,:)==i)) cycle
        if (Parent(i,1)/=0 .and. Parent(i,2)/=0 .and. Parent(i,1)==Parent(j,1) .and. &
        Parent(i,2)==Parent(j,2)) cycle  ! full sibs.  
        if (Sex(j)/=3) then
            call CalcPair(i, j, Sex(j), .FALSE., LLtmp(:,1), 1)
        else
            call CalcPair(i, j, 1, .FALSE., LLtmp(:,1), 1)
        endif
        if (Sex(i)/=3) then
            call CalcPair(j, i, Sex(i), .FALSE., LLtmp(:,2), 1)
        else
            call CalcPair(j, i, 1, .FALSE., LLtmp(:,2), 1)
        endif
        do k=1,7
            LLtmp(k,3) = MaxLL(LLtmp(k,1:2)) 
        enddo
        call BestRel(LLtmp(:,3), 1, topX, dLL)
        maybe = 1
        if (topX>3) then
            do k=1,2
                if (Parent(i,k)>0) then
                    if (ANY(Parent(Parent(i,k),:)==j))  maybe = 0  ! GP
                    if (Parent(i,k) == Parent(j,k))  maybe = 0
                endif
                if (Parent(j,k)>0) then
                    if (ANY(Parent(Parent(j,k),:)==i))  maybe = 0
                endif
            enddo
        endif
        if (maybe == 1) then
            nAmbiguous = nAmbiguous + 1
            write (201,'(2a15, 2i5, i10, a10, 2f10.2)') Id(i), ID(j), Sex(i), Sex(j), &
              AgeDiff(i,j), RelName(TopX), LLR_O(i,j), dLL
        endif
    enddo
enddo   
 close (201)

deallocate(BYRank)
deallocate(OppHomDF)
deallocate(SortBY)
 call DeAllocAll
 
end subroutine parents
!end program AssignParents

! ##############################################################################################

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################

subroutine sibships
!program Sibships
use Global
implicit none

integer :: Round, i, j, s, k, n, RX, LastR, topX
 character(len=1) :: RoundC
 character(len=2) :: RelName(8)
 character(len=4) :: DumTmp
 character(len=30), allocatable, dimension(:,:) :: DumName
 character(len=30), allocatable, dimension(:,:,:) :: GpName
double precision :: TotLL(42), dLL, tmpLL(7), LRS

RX = 1  ! no. of initial rounds, pairs-cluster-merge only
    
 call ReadData
 call ReadParents
 call PrecalcProbs
 
allocate(PairID(5*nInd, 2)) 
allocate(PairDLLR(5*nInd))
allocate(PairType(5*nInd))  ! maternal (1), paternal (2) or unknown (3)
allocate(CLL(nInd/2,2))
 CLL = 999
allocate(nS(nInd/2,2))
nS = 0
allocate(SibID(maxSibSize, nInd/2, 2))
SibID = 0
allocate(GpName(2,nInd/2,2))
allocate(DumName(nInd/2, 2))
allocate(LR_parent(nInd,3))
LR_parent = 999
allocate(LR_GP(3, nInd/2,2))
LR_GP = 999
nC = 0

!=======================

! find current FS (based on real parents)
do i=1,nInd-1
    do j=i,nInd
        if (Parent(i,1)==Parent(j,1) .and. Parent(i,1)/=0 .and. &
            Parent(i,2)==Parent(j,2) .and. Parent(i,2)/=0) then
            call MakeFS(i, j)
        endif
    enddo
enddo
 call UpdateAllProbs
! call dblepr("Sibships - Initial Total LL : ", -1, SUM(Lind), 1)
call intpr ( "Sibships ... ", -1, 0, 0)
 
LastR = 0
do Round=1, Nrounds
    TotLL(Round) = SUM(Lind)
    write(RoundC, '(i1)') Round  ! for filenames
    if(Round > 1) then
        if (ABS(TotLL(Round-1) - TotLL(Round)) < thLR) then
            LastR = 1
        endif
    endif
        
    if (Round==1 .and. nAgeClasses>1) then
        call FindPairs(.FALSE.)   ! do not use age priors yet to diff. HS/GP/FA
    else
        call FindPairs(.TRUE.)
    endif
    
!    open(unit=405, file="PairwiseLLR_"//RoundC//".txt", status="replace")
!    write(405, '(2a15, 2a10, a15)') "ID1", "ID2", "PatMat", "dLLR"
!    do i=1,nPairs
!        write(405, '(2a15, i10, f15.4)') ID(PairID(i,1:2)), PairType(i),  PairDLLR(i)
!    enddo
!    close(405)
   
    if (Round==1) then
        call Clustering(-1)
    else
        call Clustering(LastR)
    endif
 !  call intpr ( "Clusters (mat): ",16, nC(1), 1)
 !  call intpr ( "Clusters (pat): ",16, nC(2), 1)
    call UpdateAllProbs
    
    call Merging
    call UpdateAllProbs
    
    if (Round > RX .or. Round==Nrounds .or. LastR==1) then
        call GrowClusters
        call UpdateAllProbs
    endif
   
 !   open (unit=401,file="Sibships_"//RoundC//".txt",status="replace")
 !   write (401,'(2a10, a5)')  "PatMat", "Sibship", "ID"
 !   do k=1,2
 !       do s=1,nC(k)
 !           do n=1, nS(s,k)
 !               write (401,'(2i10, " ", a30)') k, s, Id(SibID(n,s,k))
 !           enddo
 !       enddo
 !   enddo
 !   close (401)
 
    if (nAgeClasses > 1 .and. (Round > RX .or. Round==Nrounds .or. LastR==1)) then     
        call SibParent  ! replace dummy parents by indivs
        call UpdateAllProbs
        
        call MoreParent  !  assign additional parents to singletons (e.g. unknown BY/sex)
        call UpdateAllProbs
        
        call SibGrandparents
        call UpdateAllProbs
    endif
    
    !=======================
    call dblepr("Round "//RoundC//", Total LL : ", -1, SUM(Lind), 1) 
    
    ! give names to dummy parents
    do k=1,2
        do s=1, nC(k)
            write(DumTmp, '(i4.4)') s
            DumName(s,k) = trim(DumPrefix(k))//DumTmp
        enddo
    enddo    
    
    do i=1,nInd
        do k=1,2
            if (Parent(i,k)>0) then
                ParentName(i,k) = Id(Parent(i,k))
            else if (Parent(i,k)<0) then
                ParentName(i,k) = DumName(-Parent(i,k),k)
            endif
        enddo
    enddo
    
    GpName = "NA"
    do k=1,2
        do s=1,nC(k)
            do n=1,2
                if (GpID(n,s,k)>0) then
                    GpName(n,s,k) = Id(GpID(n,s,k))
                else if (GpID(n,s,k)<0) then
                    GpName(n,s,k) = DumName(-GpID(n,s,k), n)
                endif
            enddo
        enddo
    enddo
        
    if (Round == nRounds .or. LastR==1) then
        call UpdateAllProbs
        call CalcParentLLR    
        
        open (unit=407, file=trim(PedigreeFileName), status="unknown")
        write (407,'(10a15)')  "ID", "Dam", "Sire", "LLR_dam", "LLR_sire", "LLR_pair"
        do i=1,nInd
            write (407,'(3a20, 3f10.2)') Id(i), ParentName(i,1:2), LR_parent(i,:)
        enddo
        do k=1,2
            do s=1,nC(k)
                write (407,'(3a20, 3f10.2)') DumName(s,k), GpName(:, s,k), LR_GP(:, s, k)
            enddo
        enddo
        close (407)
        
        ! output non-assigned likely 1st/2nd degree relatives pairs. 
        RelName = (/ "PO", "FS", "HS", "GP", "FA", "HA", "U ", "XX" /) 
        open (unit=201, file="Unassigned_relatives_sib.txt", status="unknown")
        write (201, '(2a15, 2a5, 3a10)') "ID1", "ID2", "Sex1", "Sex2", "BestRel", "LLR_R_U", "LLR_R1_R2"
        do i=1,  nInd-1
            do j=i+1,nInd
                do k=1,2
                    if (Parent(i,k)/=0 .and. Parent(j,k)/=0) cycle  
                    ! assume assigned 1st/2nd degree relatives are not futher related.
                    if (ANY(Parent(i,:)==j) .or. ANY(Parent(j,:)==i))  cycle  ! PO
                    if (Parent(i,k)/=0 .and. Parent(i,k)==Parent(j,k)) cycle
                    if (Parent(i,k) < 0) then
                        if (ANY(GpID(:,-Parent(i,k),k)==j))  cycle
                    endif
                    if (Parent(j,k) < 0) then
                        if (ANY(GpID(:,-Parent(j,k),k)==i))  cycle  
                    endif
                    if (k==2 .and. Parent(i,1)==0 .and. Parent(j,1)==0) cycle ! already done
                    call PairQS(i, j, LRS)  ! quick check
                    if (LRS < thLRrel) cycle           
                    if (AgeDiff(i,j)==999 .or. AgeDiff(i,j)>=0) then
                        call CalcPair(i, j, k, .FALSE., tmpLL, 4)  ! 4: gp -> return all
                    else
                        call CalcPair(j, i, k, .FALSE., tmpLL, 4)
                    endif
                    call BestRel(tmpLL, 3, topX, dLL)
                    if (topX>4) cycle
                    write (201,'(2a15, 2i5, i10, a10, 2f10.2)') Id(i), ID(j), Sex(i), Sex(j), &
                      AgeDiff(i,j), RelName(TopX), LRS, dLL  
                enddo
            enddo
        enddo 
        close (201) 
        
        call intpr ( "Finished. ", -1, 0, 0)   ! Round
        exit
!    else
!        open (unit=407,file="Pedigree_"//RoundC//".txt",status="replace")
!        write (407,'(3a20)')  "ID", "Dam", "Sire"
!        do i=1,nInd
!            write (407,'(3a20,3f15.4)') Id(i), ParentName(i,1:2)
!        enddo
!        do k=1,2
!            do s=1,nC(k)
!                write (407,'(3a20,3f15.4)') DumName(s,k), GpName(:, s,k)
!            enddo
!        enddo
!        close (407)       
    endif
enddo

deallocate(DumName)
deallocate(GpName)
 call DeAllocAll
    
end subroutine Sibships
!end program Sibships

! ##############################################################################################

! @@@@   SUBROUTINES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################

subroutine CalcOppHom  ! nInd x nInd matrix with no. opp. hom. loci
use Global
implicit none

integer :: i, j, l, Lboth, OH

OppHomM = -999

do i=1, nInd-1
    if (MOD(i,500)==0) then
    endif
    do j=i+1,nInd
        OH = 0
        do l=1,nSnp
            if ((Genos(l,i)==1).and.(Genos(l,j)==3)) then
                OH = OH+1
                if (OH > maxOppHom) exit
            endif                       
            if ((Genos(l,i)==3).and.(Genos(l,j)==1)) then
                OH = OH+1
                if (OH > maxOppHom) exit
            endif                       
        enddo
        OppHomM(i,j) = OH
        OppHomM(j,i) = OH
        if (OH <= maxOppHom) then
            Lboth = COUNT(Genos(:,i)/=-9 .and. Genos(:,j)/=-9)    

            if (Lboth < nSnp/4.0) then   ! more than 3/4th of markers missing
                OppHomM(i,j) = -Lboth 
                OppHomM(j,i) = -Lboth
            endif
        endif
    enddo
enddo

end subroutine CalcOppHom

! ##############################################################################################

subroutine CalcPO(A,B, LLR)  ! LLR of A as offspring from B, vs A as random sample from pop
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LLR
integer :: l
double precision :: PrL(nSnp)

LLR = 999
PrL = 0

do l=1,nSnp
    if (Genos(l,A)/=-9 .and. Genos(l,B)/=-9) then
        PrL(l) = LOG10(OKOP(Genos(l,A), Genos(l,B), l))
    else if (Genos(l,A)/=-9) then
        PrL(l) = LOG10(OHWE(Genos(l,A),l))
    endif
enddo

LLR = SUM(PrL) - Lind(A)

end subroutine CalcPO

! ########################################################################################

subroutine Parentage(BYrank)
use Global
implicit none

integer, intent(IN) :: BYrank(nInd)  ! TODO: make global
integer :: i, j, x, y, k, CandPar(20, 2), nCP(2), u, v, CurPar(2), maxnp

maxnp = 20   ! set equal to 1st dim of CandPar (used for double checking only.)

do x=1, nInd
    i = BYRank(x)
    nCP = 0
    CandPar = 0
    do y=1,nInd 
        j = BYRank(y)
        if (i==j) cycle
        if (j==parent(i,1) .or. j==parent(i,2)) cycle
        if (OppHomM(i,j) > maxOppHom .or. OppHomM(i,j)<0) cycle 
        if (LLR_O(i,j) < -thLR) cycle 
        if (AgeDiff(i,j) == 999) then
            if (ANY(Parent(j,:)==i) .or. ANY(Parent(i,:)==j)) cycle
            if (ALL(Parent(i,:)==0) .and. ALL(Parent(j,:)==0)) cycle
        else 
            if (AgeDiff(i,j) <= 0)  cycle
            if (Sex(j)==3 .and. Parent(i,1)==0 .and. Parent(i,2)==0) cycle 
        endif
        
        call CalcPOZ(i,j)  ! assigns parent as side effect
        do k=1,2
            if (Sex(j) /= 3 .and. Sex(j)/= k) cycle
            if (nCP(k) == maxnp) cycle 
            nCP(k) = nCP(k) + 1
            CandPar(nCP(k), k) = j
        enddo
    enddo
    if (Sex(i)==3 .or. BY(i)==-999) cycle
    if (nCP(1)>0 .and. nCP(2)>0 .and. (nCP(1)>1 .or. nCP(2)>1)) then  ! check if any combo was missed
        if (nCP(1)==1 .and. nCP(2)==2 .and. ANY(CandPar(1:2, 2) == CandPar(1,1))) cycle ! only 2 unique parents
        if (nCP(2)==1 .and. nCP(1)==2 .and. ANY(CandPar(1:2, 1) == CandPar(1,2))) cycle 
        CurPar = Parent(i, :)
        do u = 1, nCP(1)
            if (Sex(CandPar(u,1))==3 .or. BY(CandPar(u,1))==-999) cycle
            if (CandPar(u,1) == Parent(i,1)) cycle
            do v = 1, nCP(2)
                if (CandPar(v,2) == Parent(i,2)) cycle
                if (CandPar(u, 1) == CandPar(v, 2)) cycle ! when unknown sex
                Parent(i, 1) = CandPar(u, 1)
                call CalcPOZ(i,CandPar(v, 2))  ! removes unlikely parents
            enddo
        enddo
        
        do k=1,2
            if (Parent(i,k) == CurPar(k) .or. CurPar(k)==0) cycle
            call CalcPOZ(i, CurPar(k))  ! opportunity to restore previous parent
        enddo
    endif
enddo

end subroutine Parentage

! ########################################################################################
subroutine CalcPOZ(A, B)  ! replace a current parent of A by B? k=sex(B)
use Global
implicit none 

integer, intent(IN) :: A, B
integer :: m, CurPar(2), TopX, k
double precision :: LLA(2,7,7), TopLL, LLcp(3), LLBA(7), TopBA, dLL, LLtmp(2)
 
 CurPar = Parent(A,:)

if (AgeDiff(A,B)==999 .and. (ALL(Parent(A,:)==0) .and. ALL(Parent(B,:)==0))) then
    return  ! can't tell if A or B is parent
endif

if (Sex(B)/=3) then
    k = Sex(B) 
else if (ALL(Parent(A,:)==0) .or. ALL(Parent(A,:)/=0)) then
    return  ! can't tell if B is father or mother
else 
    do m=1,2
        if (Parent(A,m)==0) then
            k = m  ! try B as opposite-sex parent
        endif
    enddo   
endif

 call CalcLind(A)
 call CalcLind(B)

TopX = 0
LLA = 999
LLBA = 999
if (Parent(A,1)==0 .and. Parent(A,2)==0) then   
    call CalcPair(A, B, k, .FALSE., LLA(1,:,7), 1)   
    call BestRel(LLA(1,:,7), 1, TopX, dLL)
    TopLL = MaxLL(LLA(1,:,7))
    if (AgeDiff(A,B)==999) then
        if (Sex(A)/=3) then
            m = Sex(A)
        else if (Parent(B,2)==0) then
            m = 2
        else
            m = 1
        endif
        if (Parent(B,m)==0) then
            call CalcPair(B, A, m, .FALSE., LLBA, 1)
            if (Parent(B,3-m) < 0) then  ! include changes in its CLL
                call CalcU(B,m, Parent(B,3-m),3-m, LLtmp(1))
                Parent(B,m) = A
                call CalcU(B,m, Parent(B,3-m),3-m, LLtmp(2))
                Parent(B,m) = 0
                LLBA(1) = LLBA(1) + (LLtmp(2) - LLtmp(1))
            endif
            TopBA = MaxLL(LLBA)  
            if (ABS(TopBA - TopLL) < thLRrel .or. (TopBA > TopLL))  return
        endif
    endif
    if (TopX==1) then
        Parent(A,k) = B
    endif
else 
    LLcp = 0  ! need LL over all (2-4) indiv involved
    TopLL = 999
    call CalcU(CurPar(1), 1, CurPar(2), 2, LLcp(3))
    do m=1,2
        call CalcU(CurPar(3-m), 3-m, B, k, LLcp(m))
        if (curPar(3-m) < 0) then
            LLCP(m) = LLCP(m) - Lind(A)
            LLCP(3) = LLCP(3) - Lind(A)  ! never called when both parents <0
        endif
    enddo
       
    do m=1,2
        if (CurPar(m)==0) cycle  ! LLA(m,:,:) empty if only 1 CurPar
        call CalcPair(A, B, k, .FALSE., LLA(m,:,1), 1)   ! CurPar(m)=Par + B_7        
        Parent(A,m) = 0 
        if (CurPar(m) < 0)  call RemoveSib(A, -CurPar(m), m)
        
        call CalcPair(A, B, k, .FALSE., LLA(m,:,7), 1)   ! B_7
        if (CurPar(m)>0) then
            call CalcPair(A, CurPar(m), m, .FALSE., LLA(m,7,:), 1)  ! CurPar(m)_7  
            Parent(A,k) = B
            call CalcLind(A)
            call CalcPair(A, CurPar(m), m, .FALSE., LLA(m,1,:), 1)  ! B=Par, CurPar(m)_7
        else if (CurPar(m)<0) then
            if (nS(-CurPar(m),m)==1) then  ! was sibling pair
                call calcPair(A, SibID(1,-CurPar(m),m), m, .FALSE., LLA(m,7,:), 3)
            else
                call checkAdd(A, -CurPar(m), m, LLA(m,7,:), 4)  ! CurPar(m)_7
                call ReOrderAdd(LLA(m,7,:))
            endif
                
            Parent(A,k) = B
            if (Parent(A,k)<0)  call RemoveSib(A, -CurPar(k), k)
            if (nS(-CurPar(m),m)==1) then
                call calcPair(A, SibID(1,-CurPar(m),m), m, .FALSE., LLA(m,1,:), 3)
            else
                call checkAdd(A, -CurPar(m), m, LLA(m,1,:), 4)  ! B=Par, CurPar(m)_7 
                call ReOrderAdd(LLA(m,1,:))
            endif
            call CalcU(B,3-m, -CurPar(m), m, LLtmp(1))            
        endif
        if (Parent(B, m)==CurPar(m) .and.  CurPar(m)/= 0) then  ! HS implies curPar = Par
            LLA(m,3,7) = 888
        endif
        
        Parent(A,:) = CurPar  ! restore
        if (CurPar(m) < 0) call DoAdd(A, -CurPar(m), m)
        if (m/=k .and. CurPar(k)<0)  call DoAdd(A, -CurPar(k), k)
        call CalcLind(A)
        WHERE (LLA(m,2:6,1)<0) LLA(m,2:6,1) = LLA(m,2:6,1) + LLcp(3) 
        WHERE (LLA(m,2:6,7)<0) LLA(m,2:6,7) = LLA(m,2:6,7) + LLcp(3)   
        WHERE (LLA(m,7,:)<0) LLA(m,7,:) = LLA(m,7,:) + LLcp(m)     
        WHERE (LLA(m,1,:)<0) LLA(m,1,:) = LLA(m,1,:) + LLcp(m)
        if (m==3-k .and. curPar(m) < 0)  LLA(m,1,1) = LLtmp(1) + Lind(A)  ! note: CurPar(k)==0
    enddo  
    
    TopLL = MaxLL(RESHAPE(LLA(:,:,:), (/2*7*7/))) ! MAXLOC doesn't cope well with ties

    if (AgeDiff(A,B)==999) then
        if (Sex(A)/=3) then
            call CalcPair(B, A, Sex(A), .FALSE., LLBA, 1)
        else
            if (Parent(B,1)/=0 .and. Parent(B,2)/=0) return
            do m=1,2
                if (Parent(B,m)==0) then
                    call CalcPair(B, A, m, .FALSE., LLBA, 1)  ! try A as opposite-sex parent
                endif
            enddo   
        endif
        WHERE (LLBA<0) LLBA = LLBA + LLcp(3)
        TopBA = MaxLL(LLBA)
        if (ABS(TopBA - TopLL) < thLRrel .or. (TopBA > TopLL)) then
            return
        endif
    endif
    
    if (LLA(3-k,1,1)==TopLL .or. ANY(LLA(k,1,:) == TopLL)) then  ! B + CurPar(3-k) 
        Parent(A, k) = 0
        call BestRel(LLA(3-k,:,1), 1, topX, dLL)
        if (topX == 1 .and. dLL > thLRrel) then
            Parent(A, k) = B
        endif
        Parent(A, 3-k) = CurPar(3-k)
    else if ((TopLL - MaxLL(RESHAPE(LLA(:,2:7,1), (/2*6/)))) < 0.01) then  ! keep CurPar (both)
        Parent(A,:) = CurPar
!        do m=1,2
!            call BestRel(LLA(m,1,:),1,topX, TopBA)
!            if ((MaxLL(LLA(m,2:7,1)) - TopBA) > thLRrel) then
!                Parent(A,m) = CurPar(m)
!            else
 !               Parent(A,m) = 0  ! unclear
 !           endif
 !       enddo
    else if (ANY(LLA(:,1,:)==TopLL)) then  ! only B
        call BestRel(LLA(k,:,1),1,topX, TopBA)
        if ((MaxLL(LLA(k,1,:)) - TopBA) > thLRrel) then
            Parent(A,k) = B
        else
            Parent(A,k) = 0  ! unclear.
        endif
        Parent(A,3-k) = 0
    else if (ANY(LLA(k, 2:7, 2:7) == TopLL)) then  ! keep CurPar(3-k)
        Parent(A, k) = 0
    else if (ANY(LLA(3-k, 2:7, 2:7) == TopLL)) then  ! keep CurPar(k)
        Parent(A, k) = CurPar(k)
        Parent(A, 3-k) = 0
    else
        Parent(A,:) = CurPar
    endif
endif 

 call CalcLind(A)

end subroutine CalcPOZ

! ##############################################################################################

subroutine ReOrderAdd(LL)  ! reorder output from CheckAdd for compatibility with CalcPair (for POZ)
use Global
implicit none

double precision, intent(INOUT) :: LL(7)
double precision :: LLtmp(7)

LLtmp = 999
if (LL(4) - MaxLL(LL(2:3)) > -thLRrel) then  ! more likely to be parent of dummy than offspring / unclear
    LLtmp(1) = 222
else
    LLtmp(1) = MaxLL(LL(2:3))
endif
LLtmp(2:3) = LL(5:6)
LLtmp(7) = LL(7) 

LL = LLtmp

end subroutine ReOrderAdd

! ##############################################################################################

subroutine FindPairs(UseAgePrior)
use Global
use qsort_c_module
implicit none

logical, intent(IN) :: UseAgePrior
integer :: k, i, j, top, PairTypeTmp(5*nInd), PairIDtmp(5*nInd,2)
double precision :: dLL, PairLLRtmp(5*nInd), tmpLL(7), LRS
integer, allocatable, dimension(:) :: Rank
double precision, allocatable, dimension(:) :: SortB

nPairs = 0
PairID = -9
PairDLLR = 999
PairType = 0

do i=1,  nInd-1
!    if (MODULO(i,200)==0) then 
!        call intpr ( " ",1, i, 1)
!    endif
    if (Parent(i,1)/=0 .and. Parent(i,2)/=0) cycle
    do j=i+1,nInd
        do k=1,2
            if (Parent(i,k)/=0 .or. Parent(j,k)/=0) cycle
            if (k==2 .and. Parent(i,1)==0 .and. Parent(j,1)==0) cycle ! already done
            call PairQS(i, j, LRS)  ! quick check
            if (LRS < -thLR) cycle  ! true FS more likely to be HS than U 
            if (AgeDiff(i,j)==999 .or. AgeDiff(i,j)>=0) then
                call CalcPair(i, j, k, UseAgePrior, tmpLL, 3)
            else
                call CalcPair(j, i, k, UseAgePrior, tmpLL, 3)
            endif
            call BestRel(tmpLL, 3, top, dLL)
            if (top==2 .or. top==3) then
                nPairs = nPairs+1
                if (nPairs > 5*nInd) then
                    call rexit("too many pairs")
                endif
                PairID(nPairs, :) = (/ i, j /)
                PairDLLR(nPairs) = dLL
                if (k==1 .and. Parent(i,2)==0 .and. Parent(j,2)==0) then
                    pairType(nPairs) = 3  
                else
                    PairType(nPairs) = k
                endif
            endif
        enddo
    enddo
enddo
 
 ! sort by decreasing dLL
PairIDtmp = 0
PairLLRtmp = 0
allocate(Rank(nPairs))
allocate(SortB(nPairs))
Rank = (/ (i, i=1, nPairs, 1) /)
SortB = PairDLLR(1:nPairs)
 
 call QsortC(SortB, Rank(1:nPairs))
do i=1,nPairs
    PairTypeTmp(i) = PairType(Rank(nPairs-i+1))  ! decreasing order
    PairIDtmp(i,1:2) = PairID(Rank(nPairs-i+1), 1:2)  
    PairLLRtmp(i) = PairDLLR(Rank(nPairs-i+1)) 
enddo 

PairType = PairTypeTmp
PairID = PairIDtmp 
PairDLLR = PairLLRtmp
deallocate(Rank)
deallocate(SortB)

end subroutine FindPairs

! ##############################################################################################

subroutine PairQS(A, B, LR)  ! quick check, not conditioning on parents. LLR HS/U
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp, 2)

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9 .or. Genos(l,B)==-9) cycle
    PrL(l,1) = LOG10(PHS(Genos(l,A), Genos(l,B), l))
    PrL(l,2) = LOG10(PFS(Genos(l,A), Genos(l,B), l))
enddo
LR = MAXVAL(SUM(PrL, DIM=1))

end subroutine PairQS

! ##############################################################################################

subroutine CalcPair(A, B, k, InclAge, LL, focal)  ! joined LL A,B under each hypothesis
use Global
implicit none
! to consider: self, full 1st cousins, GGP  (+ double 1st cousins?)

integer, intent(IN) :: A,B,k, focal
logical, intent(IN) :: InclAge  ! include age prior y/n (always checks if age compatible)
double precision, intent(OUT) :: LL(7)  ! PO,FS,HS,GG,FAU,HAU,U
integer :: x
double precision :: LLg(7), LLtmpA(2,3), LLtmpGGP, LLCC, LRS, ALR, LLX(2)

LLg = 999
LL = 999
LLtmpGGP = 999
LRS = 999

if (AgeDiff(A, B)/=999) then
    if (AgePriorM(ABS(AgeDiff(A, B))+1, k) == 0.0) then
        LLg(2) = 777  
        LLg(3) = 777  ! not sibs
    endif
    if (AgePriorM(ABS(AgeDiff(A, B))+1, 3-k) == 0.0) then
        LLg(2) = 777  
    endif
    if (AgePriorM(ABS(AgeDiff(A, B))+1,6)==0.0) then   
        LLg(5)=777
        LLg(6)=777
    endif
    if (AgeDiff(A,B) <= 0) then  ! B younger than A
        LLg(1) = 777
        LLg(4) = 777
    else if (Sex(B) /= 3 .and. Sex(B)/=k) then
        LLg(1) = 777
    else if (AgePriorM(AgeDiff(A,B)+1, k+2) == 0.0 .and. &
        AgePriorM(AgeDiff(A,B)+1, 5) == 0.0) then
        LLg(4) = 777 ! not GP
    endif
endif
if (LLg(focal) == 777) then
    LL = LLg
    return
endif

 call CalcU(A,k,B,k, LLg(7))
 
! PO?
if (LLg(1)==999) then
    call PairPO(A, B, k, LLg(1))
endif
if (focal==2 .or. focal==3 .and. Sex(B)/=k .and. Sex(B)/=3 .and. AgeDiff(A,B)>0) then
    call PairPO(A, B, Sex(B), LLg(1))   
endif

! FS? 
if (LLg(2)==999) then
    call PairFullSib(A, B, LLg(2))  ! includes check for parent mismatch
endif

! HS?
if (LLg(3)==999) then
    if (Parent(A, 3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
        LLg(3) = 777
    else
        call PairHalfSib(A, B, k, LLg(3))  ! includes check for parent mismatch
    endif
endif

! GP?
if (LLg(4)==999) then
    call PairGP(A, B, k, focal, LLg(4))
    call PairGGP(A, B, k, LLtmpGGP)
endif

! FA/HA?
if (LLg(5)==999) then  ! FA & HA have same ageprior.
    LLtmpA = 999
    do x=1,3  ! mat, pat, FS
        call PairUA(A, B, k, x, LLtmpA(1,x))
        call PairUA(B, A, k, x, LLtmpA(2,x))
    enddo
    LLg(5) = MaxLL(LLtmpA(:,3))
    LLg(6) = MaxLL(RESHAPE(LLtmpA(:,1:2), (/2*2/) ))
endif

LLCC = 999
 call PairCC(A, B, k, LLCC) 

LL = LLg
if (AgeDiff(A,B)/=999 .and. InclAge) then
    LL(1) = LLg(1)  ! TODO: age prior?
    if (LLg(2) < 0) then
        LL(2) = LLg(2) + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, 1)) &
                + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, 2)) 
    endif
    if (LLg(3) < 0) then
        LL(3) = LLg(3) + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, k))
    endif
    do x=5,6  ! same age prior for full & half aunts/uncles
        if (LLg(x) < 0) then
            LL(x) = LLg(x) + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, 6)) 
        endif
    enddo
endif

LL(6) = MaxLL( (/LL(6), LLtmpGGP, LLCC/) )  ! use most likely 3rd degree relative

if (LLg(4) < 0 .and. AgeDiff(A,B)/=999) then
        x = 0
        if (k==1) then
            if ((AgeDiff(A,B) > 0 .and. Sex(B)==1) .or. &
            (AgeDiff(A,B) < 0 .and. Sex(A)==1)) then
                x = 3  ! mat. grandmother
            else 
                x = 5  ! mat. grandfather
            endif
        else if (k==2) then
            if ((AgeDiff(A,B) > 0 .and. Sex(B)==2) .or. &
            (AgeDiff(A,B) < 0 .and. Sex(A)==2)) then
                x = 4  ! pat. grandfather
            else 
                x = 5  ! pat. grandmother
            endif
        endif
        ALR = AgePriorM(ABS(AgeDiff(A, B))+1, x)
        if (ALR /= 0) then
        if (InclAge) then
            LL(4) = LLg(4) + LOG10(ALR)
        else
            LL(4) = LLg(4)
        endif
    else
        LL(4) = 777
    endif
endif 

LLX = 999
if (LL(2)<0 .and. (focal==2 .or. focal==3)  .and. (LL(2) - MaxLL(LL)) < thLRrel) then
    if (ALL(Parent(A,:)==0) .and. ALL(Parent(B,:)==0)) then   ! consider if HS + AU
        call PairHSHA(A, B, LLX(1))
        call PairHSHA(B, A, LLX(2))
        if (MaxLL(LLX) - LL(2) > thLRrel) then
            LL(2) = 222
        endif
    endif
    if (LL(2)/=222) then ! check if inbred FS (& boost LL)
        call PairFSHA(A, B, k, LLX(1))
        call PairFSHA(A, B, 3-k, LLX(2))
        if (MaxLL(LLX) > LL(2)) then
            LL(2) = MaxLL(LLX)
        endif
    endif
endif
end subroutine CalcPair

! ##############################################################################################

subroutine PairPO(A, B, k, LL)
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp,2), PrX(3,3,2), PrPA(3), PrB(3) 

LL = 999
if(Parent(A,k)>0) then  ! allow dummy to be replaced (need for AddFS)
    if (Parent(A,k)==B) then
        LL = 888
    else
        LL = 777
    endif
endif

if (Sex(B)/=3 .and. Sex(B)/=k) then
    LL = 777
endif
if (LL/=999) return

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9) cycle
    if (Genos(l,B)==-9) then
        PrL(l,:) = LindX(l,A)
        cycle
    endif
    call ParProb(l, Parent(A,3-k), 3-k, A,0, PrPA)
    call ParProb(l, B, k, 0,0, PrB)
    
    do x=1,3
        do y=1,3
            PrX(x,y,1) = OKA2P(Genos(l, A), x, y, l) * PrB(x) * PrPA(y)
            if (Parent(A,3-k)==0) then  ! consider close inbreeding
                PrX(x,y,2) = OKA2P(Genos(l, A), x, y, l) * PrB(x) * AKAP(y,x,l)
            endif
        enddo
    enddo
    PrL(l,1) = LOG10(SUM(PrX(:,:,1)))
    if (Parent(A,3-k)==0) then
        PrL(l,2) = LOG10(SUM(PrX(:,:,2)))
    endif
enddo

if (Parent(A,3-k)==0) then
    LL = MaxLL(SUM(PrL, DIM=1)) + Lind(B)
else
    LL = SUM(PrL(:,1)) + Lind(B)
endif

end subroutine PairPO

! ##############################################################################################

subroutine PairFullSib(A, B, LL)
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: x, y, l,k, Par(2), AncA(2,16), AncB(2,16)
double precision :: PrL(nSnp), PrXY(3,3), Px(3,2), LUX(2), LLtmp

LL = 999
Par = 0  ! joined parents of A & B
if (Parent(A,1)==Parent(B,1) .and. Parent(A,1)/=0 .and. &
    Parent(A,2)==Parent(B,2) .and. Parent(A,2)/=0) then ! already FS
    LL = 888
else 
    do k=1,2
        if (Parent(A,k) == B .or. Parent(B,k) == A) then
            LL = 777
            exit
        else if (Parent(A,k)/=Parent(B,k) .and. ((Parent(A,k)>0 .and. Parent(B,k)>0) .or. &
          (Parent(A,k)<0 .and. Parent(B,k)<0))) then
            LL = 777
            exit
        else if (Parent(A,k)/=0) then
            Par(k) = Parent(A,k)
        else
            Par(k) = Parent(B,k)
        endif        
    enddo
endif  
if (LL/=999) return

 call GetAncest(A,1,AncA)
 call GetAncest(B,1,AncB)
if (ANY(AncA == B) .or. ANY(AncB == A)) then
    LL = 777
    return
endif
   
PrL = 0 
if (Par(1) < 0 .and. Par(2)<0) then  ! call AddFS
    call CalcU(A, 1, B, 1, LUX(1))
    do k=1,2 
        if (Parent(A,k)==Par(k) .and. Parent(B,k)==0) then
            call CalcU(B, k, Par(k), k, LUX(2))
            call addFS(B, -Par(k), k, 0, k, LLtmp) 
            LL = LLtmp - LUX(2) + LUX(1)
            exit
        else if (Parent(B,k)==Par(k) .and. Parent(A,k)==0) then
            call CalcU(A, k, Par(k), k, LUX(2))
            call addFS(A, -Par(k), k, 0, k, LLtmp) 
            LL = LLtmp - LUX(2) + LUX(1)
            exit
        endif
    enddo
      
else  

    do l=1, nSnp
        do k=1,2
            if (Parent(A,k)==Parent(B,k)) then  ! doesn't matter if Par(k)==0
                call ParProb(l, Par(k), k, A, B, Px(:,k))
            else if (Parent(A,k)==Par(k)) then
                call ParProb(l, Par(k), k, A, 0, Px(:,k))
            else if (Parent(B,k)==Par(k)) then
                call ParProb(l, Par(k), k, B, 0, Px(:,k))
            else
                call ParProb(l, Par(k), k, 0, 0, Px(:,k))
            endif       
        enddo 
    
        do x=1,3
            do y=1,3
                PrXY(x,y) = Px(x,1) * Px(y,2)
                if(Genos(l,A)/=-9) then
                    PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y,l)
                endif
                if(Genos(l,B)/=-9) then
                    PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,B), x, y,l)
                endif
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXY))
    enddo
    LL = SUM(PrL)
endif    

end subroutine PairFullSib

! ##############################################################################################

subroutine PairHalfSib(A, B, k, LL)
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: x, l, Par
double precision :: PrL(nSnp), PrX(3), PrPx(3,2)

LL = 999
Par = 0  ! parent K

if (Parent(A,k)/=0) then
    if (Parent(A,k)/=Parent(B,k) .and. Parent(B,k)/=0) then
        LL = 777 ! mismatch
    else if (Parent(A,k)==Parent(B,k)) then
        LL = 888 ! already captured under H0
    else
        Par = Parent(A,k)
        if (Par>0) then
            if (AgeDiff(B, Par) <= 0) then  ! Par(k) younger than B
                LL = 777
            endif
        endif
    endif
else if (Parent(B,k)/=0) then
    Par = Parent(B,k)
    if (Par>0) then
        if (AgeDiff(A, Par) <= 0) then  ! Par(k) younger than A
            LL = 777
        endif
    endif
endif
if (LL/=999) return
 
PrL = 0
do l=1,nSnp
   if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    if (Par==Parent(A,k) .and. Par/=0) then
        call ParProb(l, Par, k, A, 0, PrX)
    else if (Par==Parent(B,k) .and. Par/=0) then
        call ParProb(l, Par, k, B, 0, PrX)
    else
        call ParProb(l, Par, k, 0, 0, PrX)    
    endif
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPx(:,1))
    call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPx(:,2))
    do x=1,3
        if (Genos(l,A)/=-9) then
            PrX(x) = PrX(x) * SUM(OKA2P(Genos(l,A), x, :, l) * PrPx(:,1))
        endif
        if (Genos(l,B)/=-9) then
            PrX(x) = PrX(x) * SUM(OKA2P(Genos(l,B), x, :, l) * PrPx(:,2))
        endif
    enddo
    PrL(l) = LOG10(SUM(PrX))
enddo
LL = SUM(PrL)

end subroutine PairHalfSib

! ##############################################################################################

subroutine pairHSHA(A, B, LL)   ! HS via k, & parent A is HS of B via 3-k
use Global
implicit none

integer, intent(IN) :: A,B
double precision, intent(OUT) :: LL
integer :: l, x, y, z
double precision :: PrL(nSnp), PrXYZ(3,3,3)

if (ANY(Parent(A,:)/=0) .or. ANY(Parent(B,:)/=0)) then
    LL = 777
    return
endif   ! else not necessary.

PrL = 0
do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    do x=1,3
        do y=1,3    
            do z=1,3
                PrXYZ(x,y,z) = AHWE(y,l) * AHWE(z,l) * AKAP(x, y, l)
                if (Genos(l,B)/=-9) then
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), y, z, l)
                endif
                if (Genos(l,A)/=-9) then
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,A), x, z, l)
                endif
            enddo
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine pairHSHA

! ##############################################################################################

subroutine pairFSHA(A, B, k, LL)   ! inbred FS: parent k offspring of parent 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3), PrY(3)

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then
    LL = 444   ! TODO (prob. not necessary)
    return
endif  

PrL = 0
do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    call ParProb(l, Parent(A,3-k), 3-k, -1,0, PrY)  ! TODO: remainder of sibship
    do x=1,3
        do y=1,3    
            PrXY(x,y) = PrY(y) * AKAP(x, y, l)
            if (Genos(l,B)/=-9) then
                PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,B), x, y, l)
            endif
            if (Genos(l,A)/=-9) then
                PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y, l)
            endif
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine pairFSHA

! ##############################################################################################

subroutine PairGP(A, B, k, focal, LL)  ! calculates LL that B is maternal(k=1) or paternal(k=2) grandparent of A
use Global
implicit none

integer, intent(IN) :: A,B,K, focal
double precision, intent(OUT) :: LL
integer :: l, x, y, curGP(2), m, z, AncB(2, 16)
double precision :: PrL(nSnp,2), PrX(3,2), PrPA(3), PrG(3), LLtmp, &
PrYZ(3,3,2), PrB(3), PrLU(nSnp), PrU(3,3,3)

LL = 999
 curGP = 0  
if (Parent(A,k)>0) then
    curGP = Parent(Parent(A,k),:) ! current grandparents of A (shortcut)
else if (Parent(A,k)<0) then
    curGP = GpID(:, -Parent(A,k), k)
endif

if (AgeDiff(A, B) <= 0) then  ! B younger than A
    LL = 777
else if (Sex(B)/=3) then
    if (curGP(Sex(B))/=0) then
        if (curGP(Sex(B))/=B) then
            LL = 777  ! conflict
        endif
    endif
    m = Sex(B)
else 
    if (curGP(1)/=0 .and. curGP(2)/=0) then
        do y=1,2
            if (curGP(y) == B) then
                LL = 888
                exit
            else
                LL = 777
            endif
        enddo
    endif
    if (curGP(1)==0) then
        m = 1  ! doesn't really matter.
    else
        m = 2 
    endif
endif

if (Parent(A,k)>0 .and. LL==999) then
    if (AgeDiff(Parent(A,k), B) <= 0) then  ! B younger than Parent(A,k)
        LL = 777 
    else
        if (Sex(B)/=3) then
            call PairPO(Parent(A,k), B, Sex(B), LLtmp)
        else
            call PairPO(Parent(A,k), B, 1, LLtmp)
        endif
        if (LLtmp > 0) then    ! not exact when Parent(A,k) has low CR
            LL = 777
        endif
    endif
endif

 call GetAncest(B, k, AncB)
if (ANY(AncB == A)) then
    LL = 777
endif

if (LL/=999) return

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9) cycle
    call ParProb(l, curGP(3-m), 3-m, 0,0, PrG)
    call ParProb(l, Parent(A,3-k), 3-k, A,0, PrPA)
    call ParProb(l, B, 0, 0,0, PrB)    

    PrX = 1
    do x=1,3  ! PA(k)
        if (Parent(A,k)>0) then
            if (Genos(l,Parent(A,k))/=-9) then
                PrX(x,:) = PrX(x,:) * AcO(x, Genos(l,Parent(A,k)), l)
            endif
        else if (Parent(A,k)<0) then
            PrX(x,:) =  PrX(x,:) * DumP(x,l, -Parent(A,k),k)  ! not sure; ideally not called. 
        endif
        
        PrYZ = 0
        do y=1,3  ! PA(3-k)
            do z=1,3  !  PrG(3-m)
                PrYZ(y,z,:) = OKA2P(Genos(l,A), x, y, l) * PrG(z)
                if (Parent(A,3-k) < 0 .and. Parent(A,k)==0) then
                    PrU(x,y,z) = PrYZ(y,z,1) * PrPA(y) * AKAP(x,z,l)
                endif
                PrYZ(y,z,1) = PrYZ(y,z,1) * PrPA(y) * SUM(AKA2P(x, :, z) * PrB)  !non-inbred
                if (Parent(A,3-k)<=0) then 
                    if (Parent(A,3-k) < 0 .and. Parent(A,k)==0) then
                        if (GpID(1, -Parent(A,3-k),3-k)==0 .and. GpID(2, -Parent(A,3-k),3-k)==0 .and. &
                        Parent(B, 3-k)/=Parent(A,3-k)) then
                            PrYZ(y,z,2) =  PrYZ(y,z,2) * OKAP(Genos(l,A),y,l) * & 
                            SUM(AKA2P(x, :, z) * XPr(1,y,l, -Parent(A,3-k), 3-k) / AKAP(y,:,l) * PrB)
                        else
                            PrYZ(y,z,2) =  PrYZ(y,z,2) * SUM(AKA2P(x, :, z) * AKAP(y,:,l) * PrB)
                        endif
                    else
                        PrYZ(y,z,2) =  PrYZ(y,z,2) * SUM(AKA2P(x, :, z) * AKAP(y,:,l) * PrB)  !inbreeding loop
                    endif
                endif
            enddo
        enddo
        do z=1,2  ! non-inbred / inbred 
            PrX(x,z) = SUM(PrYZ(:,:,z))
        enddo  
    enddo
    PrL(l,:) = LOG10(SUM(PrX, DIM=1))
    if (Parent(A,3-k) < 0 .and. Parent(A,k)==0) then
        PrLU(l) = LOG10(SUM(PrU))
    endif
enddo

if (Parent(A,3-k)>0 .or. focal==1) then
    LL = SUM(PrL(:,1)) + Lind(B)
else
    LL = MAXVAL(SUM(PrL, DIM=1)) + Lind(B)     
endif

if (Parent(A,3-k) < 0 .and. Parent(A,k)==0) then
    LL = LL - SUM(PrLU) + Lind(A)  ! reference LL: Lind(A)+Lind(B)
endif

end subroutine PairGP

! ##############################################################################################

subroutine PairGGP(A, B, k, LL)   ! NOTE: code largely identical to PairGP; merge?
! calculates LL that B is maternal(k=1) or paternal(k=2) great-grandparent of A
! todo: merge with GGP for clusters
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, AncB(2, 16)
double precision :: PrL(nSnp), PrXY(3,3), LLtmp, PrPA(3,2)

LL = 999
if (AgeDiff(A, B) <= 0) then  ! B younger than A
    LL = 777
else if (B==Parent(A,k)) then
    LL = 777
else
    call GetAncest(B, k, AncB)
    if (ANY(AncB == A)) then
        LL = 777
    else if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
        LL = 444   ! unlikely. TODO: implement below.
    endif
endif
if (LL==777) return

if (Parent(A,k)>0) then 
    if (ANY(Parent(Parent(A,k), :)/=0)) then
        LL = 444    ! should be picked up elsewere
    else
        call PairGP(Parent(A,k), B, k, 4, LLtmp)
        if (LLtmp > 0) then    
            LL = LLtmp
        endif
    endif
else if (Parent(A,k)<0) then
    if (ANY(GpID(:,-Parent(A,k),k)/=0)) LL = 444
endif
if (LL/=999) return

PrL = 0    
do l=1,nSnp
    if (Genos(l,A)==-9) cycle
    if (Genos(l,B)==-9) then
        PrL(l) = LindX(l,A)
        cycle
    endif
    call ParProb(l, Parent(A,k), k, A,0, PrPA(:,k))
    call ParProb(l, Parent(A,3-k), 3-k, A,0, PrPA(:,3-k))
    do x=1,3  
        do y=1,3    ! TODO: PrY
            PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:,l) * PrPA(:,3-k)) * PrPA(x,k) * &
              AKAP(x,y,l) * AKOP(y, Genos(l,B), l)
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL) + Lind(B)
 
end subroutine PairGGP

! ##############################################################################################

subroutine PairUA(A, B, kA, kB, LL)
! B half sib or full sib (kB=3) of parent kA of A?
use Global
implicit none

integer, intent(IN) :: A,B,kA, kB  ! kB=3 : full sibs
double precision, intent(OUT) :: LL
integer :: l, x, g, y, z, GG(2), GA(2), PB(2), PA, i, nA, r, u,j,e,Ei,v, &
  AncA(2,16), AncG(2, 2,16), AA(maxSibSize), cat(maxSibSize), &
  doneB(maxSibSize), BB(maxSibSize), nB
double precision :: PrL(nSnp), PrG(3,2), PrXYZ(3,3,3), PrPA(3), PrA(3), &
  PrPB(3), PrGA(3), PrAB(3,3,3,2), PrE(3,2), PrH(3)

if (A>0) then  
    nA = 1
    AA(1) = A
    PA = Parent(A,kA)
    if (PA<0) then
        LL = 444  ! possible but not implemented; TODO
        return
    endif
else
    nA = nS(-A, kA)
    AA(1:nA) = SibID(1:nA, -A, kA)
    PA = A
endif

if (B > 0) then
    nB = 1
    BB(1) = B
    PB = Parent(B,:)
else if (B < 0) then
    nB = nS(-B, kB)
    BB(1:nB) = SibID(1:nB, -B, kB)
    PB(kB) = B
    PB(3-kB) = 0
endif

! prep
LL = 999

 call GetAncest(A, kA, AncA)
GA = AncA(:, kA+2)

if (kB < 3) then
    if (B < 0 .and. GA(kB) == B) then
        LL = 888
    else if (B > 0 .and. GA(kB) == PB(kB) .and. PB(kB)/=0) then
        LL = 888
    endif
else if (GA(1)==PB(1) .and. GA(2)==PB(2) .and. GA(1)/=0 .and. GA(2)/=0) then  ! kB==3
    LL = 888
endif

GG = 0  ! parent of B, GP of A
AncG = 0
do x=1,2
    if (x/=kB .and. kB/=3) cycle
    if (GA(x)==0) then
        GG(x) = PB(x)
    else if (GA(x)/=PB(x) .and. PB(x)/=0) then
        LL = 777
    else
        GG(x) = GA(x)
    endif
    if (ANY(AA(1:nA)==GG(x))) then
        LL = 777
    endif
    call GetAncest(GG(x), x, AncG(x, :, :))
enddo

if (A > 0) then
    if (ANY(AncG == A)) then
        LL = 777
    endif
else if (A < 0) then
    if (ANY(AncG(:, kA, 2:16) == A)) then
        LL = 777
    endif
endif
if (B > 0) then
    if (ANY(AncG == B)) then
        LL = 777
    endif
else if (B < 0) then
    if (ANY(AncG(:, kB, 3:16) == B)) then
        LL = 777
    endif
endif
if (kB<3) then
    if (B<0) then
        if (ANY(AncA(kB,3:16)==B))  LL = 444  ! B is GGP; possible but not (yet) implemented
    else if (B>0) then
        if (ANY(AncA(:,3:16)==B))  LL = 444
    endif
endif
if (LL /= 999) return

do x=2,16
    do y=1,2
        do g=1,2
            if (AncG(g,y,x) > 0) then
                if (A > 0) then
                    if (AgeDiff(A, AncG(g,y,x)) < 0) then
                        LL = 777  ! A older than putative ancestor
                    endif 
                else if (A<0) then
                    if (ANY(AgeDiff(SibID(1:nS(-A,kA),-A,kA), AncG(g,y,x)) < 0)) then
                        LL = 777  ! A older than putative ancestor
                    endif
                endif
                if (x==2) cycle  ! parent of A doesn't need to be older than B
                if (B > 0) then
                    if (AgeDiff(B, AncG(g,y,x)) < 0) then
                        LL = 777  
                    endif 
                else if (B<0) then
                    if (ANY(AgeDiff(SibID(1:nS(-B,kB),-B,kB), AncG(g,y,x)) < 0)) then
                        LL = 777 
                    endif
                endif
            endif
        enddo
    enddo
enddo
if (LL /= 999) return

!==============================================

if (A>0 .and.  B>0) then  ! quicker.
    if (ALL(Parent(A,:)>=0) .and. ALL(Parent(B,:)>=0)) then
        PrL = 0
        do l=1, nSnp
            if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
            call ParProb(l, Parent(A,3-kA), 3-kA, A,0, PrPA)
            if (kB == 3) then
                do g=1,2
                    call ParProb(l, GG(g), g, 0,0, PrG(:,g))  ! >=0
                enddo        
            else
                call ParProb(l, GG(kB), kB, 0,0, PrG(:,kB))
                call ParProb(l, GA(3-kB), 3-kB, PA,0, PrGA)
                call ParProb(l, PB(3-kB), 3-kB, B,0, PrPB)
            endif
            
            PrXYZ = 0
            do z=1,3
                do y=1,3
                    do x=1,3
                        if (kB == 3) then
                            PrXYZ(x,y,z) = AKA2P(x,y,z) * PrG(y, 1) * PrG(z, 2)
                        else
                            PrXYZ(x,y,z) = AKA2P(x,y,z) * PrG(y, kB) * PrGA(z)   
                        endif
                        if (Genos(l,A)/=-9) then
                            PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(OKA2P(Genos(l,A), x, :, l) * PrPA)  
                        endif
                        if (Genos(l,B)/=-9) then
                            if (kB==3) then
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), y, z, l)   
                            else
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(OKA2P(Genos(l,B), y, :, l) * PrPB)
                            endif
                        endif
                        if (Parent(A,kA)>0) then
                            if (Genos(l,Parent(A,kA))/=-9) then
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * AcO(x, Genos(l,Parent(A,kA)), l)
                            endif
                        endif
                    enddo
                enddo
            enddo
            PrL(l) = LOG10(SUM(PrXYZ))
        enddo
        LL = SUM(PrL)
        return
    endif
endif

!==============================================

if (A>0 .and. PA<0) then
    nA = nS(-PA, kA)
    AA(1:nA) = SibID(1:nA, -PA, kA)
endif
if (kB==3) then 
    if (B>0 .and. PB(3-kA)<0) then
        nB = nS(-PB(3-kA), 3-kA)
        BB(1:nB) = SibID(1:nB, -PB(3-kA), 3-kA)  
    endif    
else if (kB < 3) then
    if (B>0 .and. PB(kB)<0) then
        nB = nS(-PB(kB), kB)
        BB(1:nB) = SibID(1:nB, -PB(kB), kB)
    endif
endif

 cat=0
do i = 1, nA
    if (kA/=kB .and. GG(3-kA)<0 .and. Parent(AA(i), 3-kA) == GG(3-kA)) then  !incl. kB=3
        cat(i) = 1  
    else if (kA==kB .and. Parent(AA(i), 3-kA)<0) then
        if (Parent(AA(i), 3-kA)==GA(3-kA)) then
            cat(i) = 2
        else
            do j=1, nB
                if (AA(i) == BB(j)) cycle
                if (Parent(AA(i), 3-kA) == Parent(BB(j), 3-kA)) then
                    cat(i) = 3
                endif
            enddo
        endif
    endif
enddo
if (kA/=kB .and. kB/=3) then
    do j=1,nB
        if (Parent(BB(j), 3-kB) == PA .and. A>0 .and. PA<0) then
            cat(nA+1) = 4  ! only possible if A>0 (else cat(i)=1)
        endif
    enddo
endif

PrL = 0
DoneB = 0
do l=1,nSnp
    if (A>0 .and. B>0) then
        if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    endif
    do g=1,2
        if (g/=kB .and. kB/=3) cycle
        call ParProb(l, GG(g), g, -1,0, PrG(:,g)) 
    enddo
    if (kB/=3) then  ! TODO: if(ANY(cat==2))
        if (ANY(cat==2)) then
            call ParProb(l, GA(3-kB), 3-kB, -1,0, PrGA)
        else
            call ParProb(l, GA(3-kB), 3-kB, 0,0, PrGA)
        endif
        if (nB==1) then
            call ParProb(l, PB(3-kB), 3-kB, B,0, PrPB)  ! only used if all cat=0
        endif
    endif

    ! === 
    
    if (ALL(cat==0) .and. ALL(GG >=0)) then
    
        if (A < 0) then
            PrA = Xpr(1,:,l, -A,kA)
        else if (A>0) then
            call ParProb(l, Parent(A,3-kA), 3-kA, A,0, PrPA)
            if (Genos(l,A)/=-9) then
                do x=1,3 
                    PrA(x) = SUM(OKA2P(Genos(l,A), x, :, l) * PrPA)
                enddo
            else
                PrA = 1 ! ? AHWE(:,l)  
            endif
            if (Parent(A,kA)>0) then
                if (Genos(l,Parent(A,kA))/=-9) then
                    PrA = PrA * AcO(:, Genos(l,Parent(A,kA)), l)  ! prob. x given obs.
                endif
            else if (Parent(A,kA)<0) then
                do x=1,3
                    PrH(x) = Xpr(1,x,l, -Parent(A,kA),kA) / SUM(OKA2P(Genos(l,A), x, :, l) * PrPA)
                enddo
                PrH = PrH / SUM(PrH)  ! total contribution from other sibs summing to 1
                do x=1,3 
                    PrA(x) = SUM(OKA2P(Genos(l,A),x,:, l) * PrPA * PrH(x))
                enddo
            endif
        endif
    
        do x=1,3  ! PA, kA
            do y=1,3  ! PrG, kB
                do z=1,3  ! PrGA, 3-kB / PrG, 3-kB
                    if (kB==3 .and. B>0) then  ! SA/PA FS of B; call AddFS for SB/=0
                        PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z) * PrG(y,3-kA) * PrG(z,kA)    
                        if (Genos(l,B)/=-9) then
                            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), y, z,l)
                        endif
                    else 
                        PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z) * PrGA(z)
                        if (B>0) then                 
                            if (Genos(l,B)/=-9) then
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * PrG(y, kB) * &
                                  SUM(OKA2P(Genos(l,B),y,:,l) * PrPB)
                            else
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * PrG(y, kB)
                            endif
                        else if (B<0) then
                            PrXYZ(x,y,z) = PrXYZ(x,y,z) * XPr(3,y,l, -B, kB)
                        endif
                    endif     
                enddo
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXYZ))
     
    else 
    
        PrAB = 0
        if (PA>0) then
            call ParProb(l, PA, kA, -1,0, PrPA)
        else
            PrPA = 1
        endif
        
        do y=1,3  ! PrG, kB
            do x=1,3  ! PA, kA
                do z=1,3
                    if (kB==3) then
                        PrAB(x,y,z,:) = PrPA(x) * AKA2P(x,y,z) * PrG(z, kA) * PrG(y, 3-kA)
                    else
                        PrAB(x,y,z,:) = PrPA(x) * AKA2P(x,y,z) * PrGA(z) * PrG(y, kB) 
                    endif
                enddo
            enddo   ! Dim4: 1:all, 2: non-sibs only
                
            do x=1,3
                doneB = 0
                do r=1, nA
                    if (A<0 .and. NFS(AA(r)) == 0) cycle  ! moved to its FS
                    if (cat(r)==0 .or. cat(r)==3) then
                        if (A<0 .or. nA==1 .or. (nA>1 .and. ANY(FSID(1:nFS(AA(r)),AA(r))==A))) then
                            call ParProb(l, Parent(AA(r), 3-kA), 3-kA, -1, 0, PrE(:,1))
                        else
                            call ParProb(l, Parent(AA(r), 3-kA), 3-kA, AA(r), -1, PrE(:,1))    
                        endif
                        PrE(:,2) = PrE(:,1)  ! 1: all, 2: non-sibs                 
                    else
                        PrE = 1
                    endif
                    
                    do e=1,3
                        do i=1, MAX(nFS(AA(r)), 1)
                            if (Genos(l,FSID(i, AA(r)))==-9) cycle
                            if (A<0 .and. nFS(AA(r))==0) cycle
                            if (A<0  .or. nA==1 .or. (nA>1 .and. FSID(i, AA(r))==A)) then 
                                PrE(e,1) = PrE(e,1) * OKA2P(Genos(l,FSID(i,AA(r))), x, e, l) 
                            else if (nA>1) then
                                PrE(e,:) = PrE(e,:) * OKA2P(Genos(l,FSID(i,AA(r))), x, e, l)
                            endif
                        enddo

                        if (cat(r)==3) then  ! kA==kB, Parent(AA(r), 3-kA)/=0   .or. cat(nA+1)==4
                            do j=1,nB
                                if (Parent(BB(j), 3-kB) /= Parent(AA(r), 3-kB)) cycle
                                if (B<0 .and. nFS(BB(j))==0) cycle
                                do i=1, MAX(nFS(BB(j)), 1)
                                    if (Genos(l,FSID(i, BB(j)))==-9) cycle
                                    if (B<0 .and. nFS(BB(j))==0) cycle
                                    if (B<0 .or. (nB==1 .and. FSID(1, BB(j))==B) .or. (nB>1 .and. FSID(i, BB(j))==B)) then
                                        PrE(e,1) = PrE(e,1) * OKA2P(Genos(l,FSID(i,BB(j))), y, e, l)
                                    else if (nB>1) then
                                        PrE(e,:) = PrE(e,:) * OKA2P(Genos(l,FSID(i,BB(j))), y, e, l)
                                    endif
                                enddo
                                DoneB(j) = 1
                            enddo
                        endif
                      
                        if (Parent(AA(r), 3-kA) < 0 .and. cat(r)/=1 .and. &
                          (A<0 .or. nA==1 .or. (nA>1 .and. ANY(FSID(1:nFS(AA(r)), AA(r))==A)))) then   ! not: Parent(AA(r), 3-kA)== PB(kB)
                            do g=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                                Ei = SibID(g, -Parent(AA(r), 3-kA), 3-kA)
                                if (nFS(Ei) == 0) cycle 
                                if (nA>1 .and. Parent(Ei, kA) == PA .and. PA/=0) cycle
                                if (kB<3) then
                                    if (nB>1 .and. Parent(Ei, kB)== PB(kB) .and. PB(kB)/=0) cycle
                                endif
                                if (nA>1 .and. ANY(FSID(1:nFS(Ei), Ei) == A)) cycle
                                if (nB>1 .and. ANY(FSID(1:nFS(Ei), Ei) == B)) cycle
                                call ParProb(l, Parent(Ei, kA), kA, Ei,-1, PrH)  ! Assume is not one of GPs
                                do i=1, nFS(Ei)
                                    if (A>0 .and. FSID(i, Ei)==A) cycle
                                    if (B>0 .and. FSID(i, Ei)==B) cycle
                                    if (Genos(l,FSID(i, Ei))/=-9) then
                                        PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e, l)
                                    endif
                                enddo
                                PrE(e,:) = PrE(e,:) * SUM(PrH)
                            enddo
                        endif
                    enddo
                          
                    do v=1,2
                        if (cat(r)==0 .or. cat(r)==3) then  ! incl. cat(nA+1)==4
                            PrAB(x,y,:,v) = PrAB(x,y,:,v) * SUM(PrE(:,v))
                        else if (cat(r)==1) then
                            PrAB(x,y,:,v) = PrAB(x,y,:,v) * PrE(y,v)
                        else if (cat(r)==2) then
                            PrAB(x,y,:,v) = PrAB(x,y,:,v) * PrE(:,v)
                        endif
                    enddo
                enddo  ! r
            enddo  ! x
                
            do j=1,nB
                if (B<0 .and. nFS(BB(j))==0) cycle
                if (DoneB(j)==1) cycle
                if (kB/=3) then
                    if (kA/=kB .and. A<0 .and. Parent(BB(j), 3-kB)==PA) cycle
                endif
                if (kB==3) then
                    PrE = 1
               else if (B<0 .or. nB==1 .or. (nB>1 .and. ANY(FSID(1:nFS(BB(j)), BB(j))==B))) then
                    call ParProb(l, Parent(BB(j), 3-kB), 3-kB, -1, 0, PrE(:,1))
                else 
                    call ParProb(l, Parent(BB(j), 3-kB), 3-kB, BB(j), -1, PrE(:,1))
                endif
                PrE(:,2) = PrE(:,1)
                
                do e=1,3
                    do u=1, MAX(nFS(BB(j)), 1)
                        if (Genos(l,FSID(u, BB(j)))==-9) cycle
                        if (B<0 .and. nFS(BB(j))==0) cycle
                        if (B<0 .or. nB==1 .or. (nB>1 .and. FSID(u, BB(j))==B)) then
                            PrE(e,1) = PrE(e,1) * OKA2P(Genos(l,FSID(u,BB(j))), y, e, l)
                        else if (nB>1) then
                            if (kA/=kB .and. ANY(AA(1:nA) == FSID(u, BB(j)))) cycle
                            PrE(e,:) = PrE(e,:) * OKA2P(Genos(l,FSID(u,BB(j))), y, e, l)
                        endif
                    enddo
                    
                
                    if (kB/=3) then  ! TODO  kB==3
                        if (Parent(BB(j), 3-kB) < 0 .and. (B<0 .or. nB==1 .or. &
                      (nB>1 .and. ANY(FSID(1:nFS(BB(j)), BB(j))==B)))) then                   
                            do g=1, nS(-Parent(BB(j), 3-kB), 3-kB)
                                Ei = SibID(g, -Parent(BB(j), 3-kB), 3-kB)
                                if (nFS(Ei) == 0) cycle
                                if (Parent(Ei,kB)/=0) then
                                    if (nB>1 .and. Parent(Ei, kB) == PB(kB)) cycle  
                                endif
                                if (nA>1 .and. Parent(Ei, kA)== PA .and. PA/=0) cycle
                                if (nA>1 .and. ANY(FSID(1:nFS(Ei), Ei) == A)) cycle
                                if (nB>1 .and. ANY(FSID(1:nFS(Ei), Ei) == B)) cycle
                                call ParProb(l, Parent(Ei, kB), kB, Ei,-1, PrH)  ! Assume is not one of GPs
                                do i=1, nFS(Ei)
                                    if (A>0 .and. FSID(i, Ei)==A) cycle
                                    if (B>0 .and. FSID(i, Ei)==B) cycle
                                    if (Genos(l,FSID(i, Ei))/=-9) then
                                        PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e, l)
                                    endif
                                enddo
                                PrE(e,:) = PrE(e,:) * SUM(PrH)
                            enddo
                        endif
                    else if (kB==3 .and. Parent(B, kA) < 0) then  ! if Parent(B,3-kA)<0, nB > 1                  
                        do g=1, nS(-Parent(B, kA), kA)
                            Ei = SibID(g, -Parent(B, kA), kA)
                            if (Ei==B .or. Ei==PA) cycle
                            if (nFS(Ei) == 0) cycle
                            if (Parent(Ei, 3-kA) == GG(3-kA) .and. GG(3-kA)/=0) cycle
                            call ParProb(l, Parent(Ei, 3-kA), 3-kA, Ei,-1, PrH) 
                            do i=1, nFS(Ei)
                                if (Genos(l,FSID(i, Ei))==-9) cycle
                                PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e, l)
                            enddo
                            PrE(e,:) = PrE(e,:) * SUM(PrH)
                        enddo                  
                    endif
                enddo
                
                do v=1,2
                    if (kB==3) then
                        PrAB(:,y,:,v) = PrAB(:,y,:,v) * PrE(y,v)
                    else
                        PrAB(:,y,:,v) = PrAB(:,y,:,v) * SUM(PrE(:,v))
                    endif
                enddo
            enddo
        enddo
        
        PrL(l) = LOG10(SUM(PrAB(:,:,:,1)/SUM(PrAB(:,:,:,2))))
    endif
enddo

LL = SUM(PrL)

end subroutine PairUA 

! ##############################################################################################

subroutine pairCC(A,B,k, LL)  ! 1st cousins
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, u, v
double precision :: PrL(nSnp), PrXY(3,3), PrUV, PrPA(3), PrPB(3), PrC(3,3)

LL = 999
if (Parent(A,k)/=0 .and. Parent(B,k)/=0) then
    if (Parent(A,k)==Parent(B,k)) then
        LL = 777
    endif
endif

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then  ! TODO: this is temp.
    LL = 444
endif
if (Parent(A,3-k) == Parent(B,3-k) .and. Parent(A, 3-k)/=0) then
    LL = 444
endif
if (LL/=999) return

PrL = 0
do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    call ParProb(l, Parent(A,3-k), 3-k, A,0, PrPA)
    call ParProb(l, Parent(B,3-k), 3-k, B,0, PrPB)
    
    do u=1,3  ! GG1
        do v=1,3  ! GG2
            PrUV = AHWE(u,l) * AHWE(v,l)
            do x=1,3  !PA
                do y=1,3    !PB
                    PrXY(x,y) = AKA2P(x,u,v) * AKA2P(y,u,v) * PrUV
                    if (Genos(l,A)/=-9) then
                        PrXY(x,y) = PrXY(x,y) * SUM(OKA2P(Genos(l,A), x, :, l) * PrPA)
                    endif
                    if (Genos(l,B)/=-9) then
                        PrXY(x,y) = PrXY(x,y) * SUM(OKA2P(Genos(l,B), y, :, l) * PrPB)
                    endif
                enddo
            enddo
            PrC(u,v) = SUM(PrXY)
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrC))
enddo

LL = SUM(PrL)

end subroutine pairCC

! ##############################################################################################

subroutine Clustering(LR)
use Global
implicit none

integer, intent(IN) :: LR  ! -1: first round, +1: last round
integer :: k, x, n, m, ij(2), sx(2), topX, PK, u
double precision :: LLtmp(7), dLL, LLx(7, 2)

do x=1, nPairs
    LLtmp = 999
    ij = PairID(x,:)
    PK = PairType(x)
    do k=1,2
        sx(1) = -Parent(ij(1),k)  ! interested in sibships, which have neg. numbers
        sx(2) = -Parent(ij(2),k)
        if (sx(1)==0 .and. sx(2)==0) then
            if (AgeDiff(ij(1),ij(2))==999) then
                call CalcPair(ij(1), ij(2), k, LR>0, LLx(:,1), 3)  ! HS/FS symmetrical, rest not.
                call CalcPair(ij(2), ij(1), k, LR>0, LLx(:,2), 3)
                do u=1,7
                    LLtmp(u) = MaxLL(LLx(u,:))
                enddo
            else if(AgeDiff(ij(1),ij(2))>=0) then
                call CalcPair(ij(1), ij(2), k, LR>0, LLtmp, 3)  ! only rely on ageprior in last round
            else
                call CalcPair(ij(2), ij(1), k, LR>0, LLtmp, 3)
            endif
            call BestRel(LLtmp, 2, topX, dLL)  
            if (PK==3 .and. LR<=0 .and. (topX==3 .or. (topX==2 .and. dLL<2*thLRrel))) cycle  ! unclear if MHS or PHS
            call BestRel(LLtmp, 3, topX, dLL)
            if (topX==2 .or. (topX==3 .and. (dLL > 2*thLRrel .or. LR>=0))) then
                nC(k) = nC(k)+1  ! new sibship (pair)
                if (Parent(ij(1), 3-k)/=0 .and. Parent(ij(1), 3-k)==Parent(ij(2), 3-k)) then
                    call MakeFS(ij(1), ij(2))
                endif                
                nS(nC(k), k) = 2
                SibID(1:2, nC(k), k) = ij
                Parent(ij(1), k) = -nC(k)
                Parent(ij(2), k) = -nC(k)
                if (Parent(ij(1), 3-k) /=0 .and. Parent(ij(1), 3-k)==Parent(ij(2), 3-k)) then
                    call MakeFS(ij(1), ij(2))
                endif
                call CalcCLL(nC(k), k)
                do n=1, 2
                    u = SibID(n, nC(k), k)
                    if (Parent(u,3-k) < 0) then
                        call CalcCLL(-Parent(u,3-k), 3-k)
                    endif                
                    call CalcLind(SibID(n, nC(k), k))
                enddo
                call CalcCLL(nC(k), k)
            endif
        else if (sx(1)>0 .and. sx(2)>0 .and. sx(1) /= sx(2)) then   ! else do nothing
            call CheckMerge(sx(1), sx(2), k,k, LLtmp, 1)
            call BestRel(LLtmp, 3,topX, dLL)              
            if (topX==1 .and. dLL > thLRrel * MAX(nS(sx(1),k), nS(sx(2),k))) then
                call DoMerge(sx(1), sx(2), k)
            endif
        else
            do m=1,2
                if (sx(m)>0 .and. sx(3-m)==0) then
                    call CheckAdd(ij(3-m), sx(m), k, LLtmp, 3)
                    call BestRel(LLtmp, 3,topX, dLL)
                    if ((topX==2 .or. topX==3) .and. dLL > nS(sx(m),k)*thLRrel) then
                       call DoAdd(ij(3-m), sx(m), k)
                    endif
                endif
            enddo
        endif
    enddo
enddo

end subroutine Clustering

! ##############################################################################################

subroutine Merging  ! check if any of the existing clusters can be merged
use Global
implicit none

integer :: k, s, r, n, topX, xr
double precision :: LLtmp(7,2), LLm(7), dLL

do k=1,2
    do s=1,nC(k)-1
        r = s
        do xr=s+1, nC(k)
            r = r + 1
            if (r > nC(k)) exit   ! possible due to merged sibships
            topX = 0
            LLtmp = 999
            call CheckMerge(s, r, k, k, LLtmp(:,1), 1)           
            do n = 1,7
                LLm(n) = MaxLL(LLtmp(n,:))
            enddo
            call BestRel(LLm, 1,topX, dLL)
            if (topX==1 .and. dLL > thLRrel * MAX(nS(s,k), nS(r,k))) then
                call DoMerge(s, r, k)
                r = r-1  ! otherwise a cluster is skipped
            endif
        enddo
    enddo
enddo

end subroutine Merging

! ##############################################################################################

subroutine GrowClusters
! for each individual, check if they can be added to any sibship clusters
use Global
implicit none

integer :: k, s, i, nMaybe, TypeM(200), ClM(200), topX, x, y
double precision :: LLtmp(7), dLL, dLLM(200)
    
do i=1, nInd
    nMaybe = 0
    dLLM = 999
    do k=1,2
        if (Parent(i,k)/=0) cycle  
        do s=1,nC(k)  ! sibships
            if (GpID(1,s,k)==i .or. GpID(2,s,k)==i) cycle          
            LLtmp = 999
            call CheckAdd(i, s, k, LLtmp, 3) 
            call BestRel(LLtmp, 3,topX, dLL)
            if ((topX==2 .or. topX==3) .and. dLL > thLRrel*nS(s,k)) then 
                nMaybe = nMaybe+1
                Typem(nMaybe) = k
                Clm(nMaybe) = s 
                dLLM(nMaybe) = dLL
            endif
        enddo
    enddo
    
    x = 0
    if (nMaybe==1) then
        x = 1
    else if (nMaybe==2) then
        if (TypeM(1)/=TypeM(2)) then
            x = MAXLOC(dLLM(1:2), DIM=1)  ! assume not identical dLL's
        else if (TypeM(1)==TypeM(2)) then
!            print *, "2 maybe sibships: ", ID(i), TypeM(1:nMaybe), ", ", Clm(1:nMaybe)
!            call CheckMerge(Clm(1), Clm(2), TypeM(1), TypeM(2), LLtmp, 0)
!            write(*,'("merge? ", 7f8.1)') LLtmp
!            print *, ""
        endif
    else  ! TODO: allow 1 of one type & =MAX, >1 of other
        x = 0
!        if (nMaybe > 1) then
!            print *, ">1 maybe sibship: ", ID(i), TypeM(1:nMaybe), ", ", Clm(1:nMaybe)
!        endif
    endif
    
    if (x /=0) then
        call DoAdd(i, Clm(x), TypeM(x))
        if (nMaybe==2 .and. TypeM(1)/=TypeM(2)) then
            y = MINLOC(dLLM(1:2), DIM=1)
            call CheckAdd(i, ClM(y), TypeM(y), LLtmp, 3)
            call BestRel(LLtmp, 3,topX, dLL)
            if ((topX==2 .or. topX==3) .and. dLL > thLRrel * nS(ClM(y), TypeM(y))) then
                 call DoAdd(i, CLM(y), TypeM(y))
            endif
        endif
    endif
enddo
 
end subroutine GrowClusters

! ##############################################################################################

subroutine SibParent  ! for each sibship, check if a parent can be assigned
use Global
implicit none

integer :: k, s, xs, i,j, n,maybe, topX, CurNumC(2), CandPar(50), nCandPar, &
  OH, Par, MaybeOpp
double precision :: LLtmp(7), dLL

 CurNumC = nC
maxOppHom = MaxMismatch - FLOOR(-nSNP * Er)   ! round up to nearest integer  
 
do k=1,2
    s = 0
    do xs=1, CurNumC(k)
        nCandPar = 0
        s = s+1
        if (s > nC(k)) exit   ! possible due to dissolved sibships after parent assigned

        do i=1,nInd
            if (Sex(i)/=k .and. Sex(i)/=3) cycle
            if (Parent(i,k)==-s) cycle
            if (GpID(1,s,k)==i .or. GpID(2,s,k)==i) cycle
            maybe=1
            do n=1,nS(s,k)
                j=SibID(n,s,k)
                if (AgeDiff(i,j) > 0 .and. AgeDiff(i,j)/= 999) then  ! candidate i younger than sib j 
                    maybe = 0
                    exit
                endif
                OH = COUNT(((Genos(:,i)==1).and.(Genos(:,j)==3)) .or. &
                (Genos(:,i)==1).and.(Genos(:,j)==3))
                if (OH > maxOppHom) then
                    maybe = 0
                    exit
                endif   
            enddo
           if (maybe==1) then               
                LLtmp = 999
                call CheckAdd(i, s, k, LLtmp, 1)
                call BestRel(LLtmp, 1, topX, dLL)
                if (topX/=1)  maybe = 0
            endif
            if (maybe==1 .and. Sex(i)==3) then  ! check if parent of opposite sex instead
                MaybeOpp = 1
                call getFSpar(s, k, Par)
                if (Par > 0)  MaybeOpp = 0
                if (Par==0 .and. ANY(Parent(SibID(1:nS(s,k), s,k),3-k)>0))  MaybeOpp = 0
                if (Par/=0 .and. Parent(i, 3-k) == Par)  MaybeOpp = 0  ! already are HS via opp parent   
                if (MaybeOpp == 1) then
                    if (Par < 0) then  ! may have more/fewer sibs
                        call CheckAdd(i, -Par, 3-k, LLtmp, 1)
                        call BestRel(LLtmp, 1, topX, dLL)
                        if (topX==1)  maybe = 0
                    else if (Par == 0) then
                        call PairPO(SibID(1, s, k), i, 3-k, LLtmp(1))   ! prob. not necessary todo with all sibs
                        call CalcU(SibID(1, s, k), k, i, 3-k, LLtmp(2))
                        if (LLtmp(1) - LLtmp(2) > thLR)  maybe = 0
                    endif
                endif
            endif
            if (maybe==1) then               
                nCandPar = nCandPar + 1
                CandPar(nCandPar) = i
            endif
        enddo
        
        if (nCandPar == 1) then  
            i = CandPar(1)
!            print *, "new sibship parent: ", k, s, ID(i), " - ", ID(SibID(1:nS(s,k), s, k))
            do n=1,nS(s,k)  ! assign parent i to all sibs in S
                Parent(SibID(n, s, k), k) = i
            enddo   
            if (Sex(i)==3)  Sex(i) = k
            do n=1,nS(s,k)
                call CalcLind(SibID(n, s, k))
            enddo
            call DoMerge(0, s, k)  !removes cluster s, shifts all subsequent ones, fixes GPs
            s = s-1  ! otherwise a cluster is skipped
!        else if (nCandPar > 1) then
!            print *, "> 1 cand parent ", k, s, ID(CandPar(1:nCandPar))
        endif
    enddo ! s
enddo ! k
    
end subroutine SibParent

! ##############################################################################################

subroutine MoreParent  ! for each individual, check if a parent can be assigned now.
use Global
implicit none

integer :: i, j, l, OH, Lboth
double precision :: LRP

do i=1, nInd
    if (ALL(Parent(i,:)/=0)) cycle  ! OPTIONAL: re-check sibship members individually.
    do j=1, nInd  ! candidate parent.
        if (i==j) cycle
        if (ANY(Parent(j,:)==i) .or. ANY(Parent(i,:)==j)) cycle
        if (AgeDiff(i,j) <= 0)  cycle
        if (Sex(j) == 3) then
            if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle 
        else
            if (Parent(i, Sex(j)) /= 0) cycle 
        endif
        OH = 0
        do l=1,nSnp
            if ((Genos(l,i)==1).and.(Genos(l,j)==3)) then
                OH = OH+1
                if (OH > maxOppHom) exit
            endif                       
            if ((Genos(l,i)==3).and.(Genos(l,j)==1)) then
                OH = OH+1
                if (OH > maxOppHom) exit
            endif                       
        enddo
        if (OH > maxOppHom) cycle  
        if (OH <= maxOppHom) then
            Lboth = COUNT(Genos(:,i)/=-9 .and. Genos(:,j)/=-9)    
            if (Lboth < nSnp/4.0)  cycle   ! more than 3/4th of markers missing
        endif        
        
        call CalcPO(i, j, LRP)
        if (LRP < -thLR) cycle
        call CalcPOZ(i,j)  ! assigns parent as side effect
    enddo
enddo

end subroutine MoreParent

! ##############################################################################################

subroutine SibGrandparents 
! for each sibship, find the most likely male and/or female grandparent
use Global
implicit none

integer :: k, s, i,j, n,maybeGP, r, m, par
double precision :: LRG 

do k=1,2
    do s=1, nC(k) 
        call CalcCLL(s,k)
        if (ALL(Parent(SibID(1:nS(s,k), s, k), 3-k) < 0)) then
            call getFSpar(s, k, par)
            if (par < 0) then
                if (nS(-par, 3-k) == nS(s,k))  cycle    ! all sibs are FS, cannot tell if mat or pat GP
            endif          
        endif
        do i=1,nInd
            maybeGP=1
            if (Parent(i,k)==-s) cycle
            if (GpID(1,s,k)==i .or. GpID(2,s,k)==i) cycle
            if (Sex(i)==3 .and. ALL(GpID(:,s,k)==0)) cycle ! cannot tell if GM or GF
            do n=1,nS(s,k)
                j=SibID(n,s,k)
                if (AgeDiff(i,j)== 999) cycle    ! BY(i) - BY(j)
                if (AgeDiff(j,i) < 0) then  ! candidate i younger than sib j 
                    maybeGP = 0
                    exit
                else 
                    if (Sex(i)==k) then
                        if (AgePriorM(AgeDiff(j,i)+1,k+2)==0.0) then
                            maybeGP = 0
                            exit
                        endif
                    else
                        if (AgePriorM(AgeDiff(j,i)+1, 5)==0.0) then
                            maybeGP = 0
                            exit
                        endif                 
                    endif
                endif 
            enddo
            if (maybeGP==0) cycle
            if(GPID(1,s,k)/=0 .and. GpID(1,s,k)==Parent(i,1) .and. &  ! FS of dummy
                GPID(2,s,k)/=0 .and. GpID(2,s,k)==Parent(i,2)) cycle       
            call Qadd(i, 0, s, k, LRG)  ! quick check
            if (LRG < -thLR) cycle
            call CalcGPz(s, k, i, 0, 0)  ! assigns GP as side-effect
        enddo
        
        do m=1,2
            do r=1, nC(m)
                if (m==k .and. s==r) cycle
                if (GPID(m,s,k) == -r) cycle  ! already assigned.
                if (GPID(k,r,m) == -s) cycle
                call Qadd(-s, k, r, m,  LRG)  ! TODO: age checks? (CalcAgeLRCAU)
                if (LRG < -thLR) cycle                               
                call CalcGPz(s, k, 0, r, m)  
            enddo
        enddo
        
    enddo
enddo
                        
end subroutine SibGrandparents

! ##############################################################################################

subroutine CalcGPz(SA, kA, B, SB, kB) ! B or SB grandparent of sibship SA?
use Global
implicit none

integer, intent(IN) :: SA, kA, B, SB, kB
integer :: m,n, CurGP(2), topX, x, mid(5), BB, y, CY(4), kY(4), v
double precision :: LLA(2,7,7), dLL, TopLL, LLcp(3,2), curCLL, LLU(3), LLtmp(3)
logical :: ConPar(4,4)

 curGP = GpID(:,SA,kA)
 dLL = 999
 TopLL = 999
 
 call CalcCLL(SA,kA)
 curCLL = CLL(SA,kA)
 
 if (SB /=0) then
    n = kB
    BB = -SB
 else if (B/=0) then
    BB = B
    if (Sex(B)/=3) then
        n = Sex(B)
    else if (GpID(1,SA,kA)==0) then
        n = 1
    else !if (GpID(2,SA,kA)==0) then
        n = 2  ! TODO: check in both configs (as in POZ)
    endif
endif

LLA = 999
if (curGP(1)==0 .and. curGP(2)==0) then
    if (B /=0) then
        call checkAdd(B, SA, kA, LLA(1,:,7), 4)
    else if (SB /= 0) then
        call checkMerge(SA, SB, kA, kB, LLA(1,:,7), 4) 
    endif
    call BestRel(LLA(1,:,7), 4,topX, dLL)
    if (topX==4) then   ! .and. dLL>thLRrel*nS(SA,kA)
        GpID(n,SA,kA) = BB
    endif       
else
    ConPar = .FALSE.
    GpID(:,SA,kA) = 0
    call CalcCLL(SA,kA)
    if (ANY(curGP<0) .or. BB<0) then
        do m=1,2
            if (curGP(m)==0) cycle
            call Connected(CurGP(m), m, -SA, kA, ConPar(4,m))
            call Connected(CurGP(m), m, BB, n, ConPar(3,m)) 
            call Connected(curGP(3-m),3-m,CurGP(m), m, ConPar(2,1))
        enddo
        call Connected(BB,n, -SA, kA, ConPar(4,3))
    endif
    
    LLU = 0
    do m=1,2
        if (curGP(m)==0) cycle
        call CalcU(curGP(m),m, -SA,kA, LLtmp(1))
        LLU(m) = LLtmp(1) - CLL(SA,kA)
    enddo  
    call CalcU(BB,n, -SA,kA, LLtmp(2)) 
    LLU(3) = LLtmp(2) - CLL(SA,kA)
    do m=1,2
        LLcp(m,m) = LLU(3-m) + LLU(3) 
    enddo
    LLcp(3,:) = LLU(1) + LLU(2)  ! B
    
    if (ANY(ConPar)) then    ! calculate 'remainder' LL for each GP
        CY = (/ CurGP(1), CurGP(2), BB, -SA /)
        kY = (/ 1, 2, n, kA /)
        do m=1,2        
            do y=1,3  ! focal: CurGP, B 
                if (y/=m .and. y/=3) cycle           
                if (y==1) then
                    call CalcU(CY(2),kY(2), CY(3),kY(3), LLCP(1,m))
                else if (y==2) then
                    call CalcU(CY(1),kY(1), CY(3),kY(3), LLCP(2,m))    
                else if (y==3) then
                    call CalcU(CY(1),kY(1), CY(2),kY(2), LLCP(3,m))
                endif
                do x=1,3
                    if (ConPar(4,x) .and. x/=y) then ! focal relationship: SA, CY(y)
                        call CalcU(CY(x),kY(x), CY(4),kY(4), LLtmp(1))
                        call CalcU(CY(x),kY(x), 0,0, LLtmp(2))
                        call CalcU(CY(4),kY(4), 0,0, LLtmp(3))
                        LLCP(y,m) = LLCP(y,m) + (LLtmp(1)-LLtmp(2)-LLtmp(3))
                    endif
                    do v=1,2
                        if (ConPar(x,v) .and. (x==y .or. v==y)) then
                            call CalcU(CY(x),kY(x), curGP(v),v, LLtmp(1))
                            call CalcU(CY(x),kY(x), 0,0, LLtmp(2))
                            call CalcU(curGP(v),v, 0,0, LLtmp(3))
                            LLCP(y,m) = LLCP(y,m) + (LLtmp(1)-LLtmp(2)-LLtmp(3))
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif

    mid = (/1,2,3,5,6/)
    GpID(:,SA,kA) = CurGP 
    call CalcCLL(SA,kA)
    dLL = 0
    LLtmp = 999
    do m=1,2 ! sex currently assigned GP
        if (CurGP(m)==0) cycle
        
        if (B /=0) then
            call checkAdd(B, SA, kA, LLA(m,:,4), 4)   ! CurGP(m)=GP + A_7 
            if (Sex(B)==3)  dLL = LLA(m,4,4)
        else if (SB /= 0) then
            call checkMerge(SA, SB, kA, kB, LLA(m,:,4), 4) 
        endif
       
        GpID(m,SA,kA) = 0
        call CalcCLL(SA,kA)
        if (B /= 0) then
            call checkAdd(B, SA, kA, LLA(m,:,7), 4)   ! A_7
            if (curGP(m)/=0 .and. Parent(B,m)==curGP(m)) then ! HS w. sib = curGP being GP
                LLA(m, 6, 7) = 888
            endif   ! TODO: double check; is this really a problem? 
            if (curGP(m)>0) then
                if (Parent(CurGP(m),n)==B)  LLA(m, 6, 7) = 888   
            else if (curGP(m) < 0) then
                if (GpID(n,-CurGP(m),m)==B)  LLA(m, 6, 7) = 888
            endif
        else if (SB /= 0) then
            call checkMerge(SA, SB, kA, kB, LLA(m,:,7), 4) 
            if (curGP(m)/=0 .and. (GpID(m,SB,kB)==curGP(m) .or. &  ! parentHFS w. sib = curGP being GP
              ANY(SibID(1:nS(SB,kB),SB,kB)==curGP(m)))) then  ! SB being GGP = curGP being GP
                LLA(m, 6, 7) = 888   ! TODO: consider other 3rd degree rel only 
            endif
        endif

        if (curGP(m) > 0) then
            call checkAdd(CurGP(m), SA,kA, LLA(m,7,:), 4)  ! CurGP(m)_7
        else if (curGP(m) < 0) then
            call checkMerge(SA, -CurGP(m), kA, m, LLA(m,7,:), 4)
            if (m/=n .and. ANY(Parent(SibID(1:nS(-curGP(m),m),-curGP(m),m),3-m)==BB)) then
                call PairUA(-SA, CurGP(m), kA, m, LLA(m,7,4))   ! not counting SA FS with a sib
            endif  ! not counting SA FS with a sib
        endif
        
        if (LLA(m,4,4)<0 .or. LLA(m,4,7)<0) then  ! Else not possible
            GpID(n,SA,kA) = BB
            call CalcCLL(SA,kA)
            if (curGP(m) > 0) then
                call checkAdd(CurGP(m),SA,kA, LLA(m,4,:), 4)  ! B=GP + CurGP(m)_7
            else if (curGP(m) < 0) then
                call checkMerge(SA, -CurGP(m), kA, m, LLA(m,4,:), 4)
            endif
        endif

        WHERE (LLA(m,mid,4)<0) LLA(m,mid,4) = LLA(m,mid,4) + LLcp(3,m)
        WHERE (LLA(m,mid,7)<0) LLA(m,mid,7) = LLA(m,mid,7) + LLcp(3,m)
        WHERE (LLA(m,4,:)<0) LLA(m,4,:) = LLA(m,4,:) + LLcp(m,m)
        WHERE (LLA(m,7,:)<0) LLA(m,7,:) = LLA(m,7,:) + LLcp(m,m)
        
        if (n==m) then
            if (B>0) then
                if (Sex(B) == 3) then
                    LLA(m,4,4) = dLL  ! from checkAdd(B, SA,..)
                else
                    LLA(m,4,4) = 777
                endif
            else   
                LLA(m,4,4) = 777 ! cannot have 2 same-sex GPs
            endif
        endif
        
        GpID(:,SA,kA) = CurGP  ! restore
        call CalcCLL(SA,kA)  
    enddo    
    TopLL = MaxLL(RESHAPE(LLA(:,:,:), (/2*7*7/)))

    if (LLA(3-n,4,4)==TopLL .or. ANY(LLA(n,4,:)==TopLL)) then
        GpID(n,SA,kA) = 0
        call BestRel(LLA(3-n,:,4), 4, topX, dLL)  
        if (topX==4 .and. dLL>thLRrel) then  ! *nS(SA,kA)
            GpID(n,SA,kA) = BB
        endif
        GPID(3-n, SA, kA) = CurGP(3-n)
    else if ((TopLL - MaxLL(RESHAPE(LLA(:,:,4), (/2*7/)))) < thLRrel) then ! keep CurGP (both)
        GpID(:,SA,kA) = curGP
    else if (ANY(LLA(:,4,:)==TopLL)) then  ! only B
        GpID(n,   SA, kA) = BB
        GpID(3-n, SA, kA) = 0
    else if (ANY(LLA(n, :, :) == TopLL)) then  ! keep CurGP(3-n)
        GpID(n,   SA, kA) = 0
        GpID(3-n, SA, kA) = curGP(3-n)
    else if (ANY(LLA(3-n, :, :) == TopLL)) then  ! keep CurGP(n)
        GpID(n,   SA, kA) = curGP(n)
        GpID(3-n, SA, kA) = 0
    endif
endif

 call CalcCLL(SA, kA)
 do n=1, nS(SA,kA)
    if (Parent(SibID(n,sA,kA),3-kA) < 0) then
        call CalcCLL(-Parent(SibID(n,sA,kA),3-kA), 3-kA)
    endif         
    call calcLind(SibID(n, SA, kA))
 enddo
 
end subroutine CalcGPz

! ##############################################################################################

subroutine Qadd(A, kA, SB, kB, LR)
use Global
implicit none

integer, intent(IN) :: A, kA, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrL(nSnp), PrX(3), PrXY(3,3)

PrL = 0
do l=1,nSnp
    if (A>0) then
        if (Genos(l,A)==-9) cycle
        do x=1,3
            PrX(x) = OKAP(Genos(l,A), x, l) * DumP(x,l,SB,kB) / AHWE(x,l)
        enddo 
        PrL(l) = LOG10(SUM(PrX))
    else
        do x=1,3
            PrX(x) = XPr(1,x,l,-A,kA) * AHWE(x,l)
            do y=1,3
                PrXY(x,y) = XPr(1,x,l,-A,kA) * AKAP(x,y,l) * DumP(y,l,SB,kB)
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXY)) - LOG10(SUM(PrX))
    endif
enddo
LR = SUM(PrL)

end subroutine Qadd

! ##############################################################################################

subroutine CheckAdd(A, SB, k, LL, focal)
use Global
implicit none

integer, intent(IN) :: A, SB, k, focal
double precision, intent(OUT) :: LL(7)
double precision :: LLg(7),  LLtmp(2,2), ALRtmp(7), LLz(3), LRHS, ALRx(2,2), LUX(3), LLM(3)
integer :: x, y, Par, MaybeOpp, i, ParTmp(2), npt
 
LL = 999
LLg = 999
LLz = 999
LRHS = 999
ALRtmp = 0

! ensure LL up to date
 call CalcCLL(SB, k)
 call CalcLind(A)
 if (Parent(A,3-k)<0) then
    call CalcCLL(-Parent(A,3-k), 3-k)
endif 
! quick check
 call Qadd(A, 0, SB, k, LRHS)  ! 2nd degree relatives vs unrelated
 if (LRHS < -nS(SB,k)*thLR .and. focal/=4) return   ! need all LL's for GPZ
if (focal==1) then
    if (Sex(A)/=3 .and. Sex(A)/=k) then
        LLg(1) = 777
        return
    endif
endif

 call CalcU(A,k, -SB, k, LLg(7))   ! unrelated
LL(7) = LLg(7)

 !=======
if (LLg(1)/=777) then
    call AddParent(A, SB, k, LLg(1))  ! A parent of SB
endif
 LL(1) = LLg(1)
 call AddFS(A, SB, k,0,k, LLg(2))
 call AddSib(A, SB, k, LLg(3))
 call CalcAgeLR(SB, k, A, 0, 0,0, ALRtmp(3))
 ALRtmp(2) = ALRtmp(3)
do x=2,3
    if (LLg(x) < 0 .and. ALRtmp(x) /= 777) then
        LL(x) = LLg(x) + ALRtmp(x)
    else
        LL(x) = 777
    endif
enddo

 call AddGP(A, SB, k, LLg(4))
if (LLg(4) < 0) then 
    call CalcAgeLR(SB, k, 0, A, 0,0, ALRtmp(4))
    if (ALRtmp(4) /= 777) then
        LL(4) = LLg(4) + ALRtmp(4)
    else
        LL(4) = 777
    endif
else
    LL(4) = LLg(4)
endif

! FAU
 call pairUA(-SB, A, k, 3, LLg(5))    ! A FS with SB?
if (Parent(A, 3-k)==0) then  
    call getFSpar(SB, k, Par)
    LLtmp = 999
    if (Par /= 0) then
        call pairUA(A, -SB, k, k, LLtmp(1,1))
        LLg(5) = MaxLL((/ LLg(5), LLtmp(1,1) /))
    endif
endif

! HAU
LLtmp = 999
do x=1,2
    call pairUA(A, -SB, x, k, LLtmp(1,x))  
    if (Parent(A,x)>0) then
        call CalcAgeLR(SB, k, Parent(A,x), 0, 0,0, ALRx(1,x))
    else
        call CalcAgeLRCAU(A, x, -SB, k, ALRx(1,x))
    endif
    call pairUA(-SB, A, k, x, LLtmp(2,x))
    call CalcAgeLRCAU(-SB, k, A, x, ALRx(2,x))
    do y=1,2
        if (LLtmp(y,x) < 0 .and. ALRx(y,x) /= 777) then
            LLtmp(y,x) = LLtmp(y,x) + ALRx(y,x)
        else
            LLtmp(y,x) = 777
        endif
    enddo
enddo

if (LLg(5)<0 .and. ALRx(2,1)/=777 .and. ALRx(2,2)/=777) then
    LL(5) = LLg(5) + ALRx(2,1) + ALRx(2,2)  ! A FS of SB
else
    LL(5) = 777
endif
  LL(6) = MaxLL(RESHAPE(LLtmp, (/2*2/) ))

if ((LL(focal)<0 .and. LL(focal)>LL(7)) .or. focal==4) then
    call AddGGP(A, SB, k, LLz(1))
    if (LLz(1) < 0 .and. ALRtmp(4)/=777 .and. LLg(4)<0) then  ! TODO not perfect, but better comparison
        LLz(1) = LLz(1) + ALRtmp(4)
    endif
    call ParentHFS(A, 0,1, SB, k,3, LLz(2))
    call ParentHFS(A, 0,2, SB, k,3, LLz(3))
    LL(6) = MaxLL((/LL(6), LLz/))
endif
            
LLM = 999    ! TODO: this piece slows it down quite a bit; do more checks
if ((MaxLL(LL)==LL(3) .or. MaxLL(LL)==LL(2)) .and. (focal==2 .or. focal==3) .and. &
  Parent(A,3-k)==0) then  ! check if not HS via opp parent only
    MaybeOpp = 1
    call getFSpar(SB, k, Par)
    if (Par > 0)  MaybeOpp = 0
    if (Par==0) then
        if (ANY(Parent(SibID(1:nS(SB,k), SB,k),3-k)>0)) then
            MaybeOpp = 0
        else   ! check if opp. parent possibly to be merged
            npt = 0
            ParTmp = 0
            do i=1, nS(SB,k)
                if (Parent(SibID(i,SB,k), 3-k)<0 .and. &
                  .not. ANY(ParTmp == Parent(SibID(i,SB,k), 3-k))) then
                    npt = npt + 1
                    if (npt > 2) then
                        MaybeOpp = 0
                        exit
                    else
                        ParTmp(npt) = Parent(SibID(i,SB,k), 3-k)
                    endif
                endif
            enddo
            if (MaybeOpp == 1 .and. npt==2) then
                call CalcU(ParTmp(1), 3-k, ParTmp(2), 3-k, LLM(1))
                call MergeSibs(-ParTmp(1), -ParTmp(2), 3-k, LLM(2))
                if (LLM(1) - LLM(2) > thLR)  MaybeOpp = 0   
            endif
        endif
    endif
    if (Par/=0 .and. Parent(A, 3-k) == Par)  MaybeOpp = 0  ! already are HS via opp parent  
    if (MaybeOpp == 1 .and. Par < 0) then
        call Qadd(A, 0, -Par, 3-k, LLM(1))  ! 2nd degree relatives vs unrelated    
        if (LLM(1) < 0)  MaybeOpp = 0
    endif
    if (MaybeOpp == 1) then
        LLM = 999
        if (Par < 0) then  ! may have more/fewer sibs
            call AddFS(A, -Par, 3-k,0,3-k, LLM(1))
            call AddSib(A, -Par, 3-k, LLM(2))
            call CalcU(A, 3-k, Par, 3-k, LLM(3))
        else if (Par == 0) then
            call PairFullSib(A, SibID(1, SB, k), LLM(1))  ! prob. not necessary todo with all sibs
            call PairHalfSib(A, SibID(1, SB, k), 3-k, LLM(2))
            call CalcU(A, 3-k, SibID(1, SB, k), 3-k, LLM(3))
        endif
        if (LLM(1) < 0 .and. LLM(1) - LLM(2) > thLRrel) then
            LL(2) = LL(2)   ! proceed.
        else if (LLM(2) < 0 .and. LLM(2) - LLM(3) > thLR) then
            LL(2:3) = 222  ! more likely to be added to opp. parent, or unclear
        endif
    endif
endif

end subroutine CheckAdd

! ##############################################################################################

subroutine Qmerge(SA, SB, k, LR)
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrL(nSnp), PrX(3), PrXY(3,3)

PrL = 0
do l=1,nSnp
    do x=1,3
        PrX(x) = XPr(1,x,l,SA,k)*XPr(1,x,l,SB,k)* AHWE(x,l)
        do y=1,3
            PrXY(x,y) = XPr(1,x,l,SA,k)*XPr(1,y,l,SB,k)* AHWE(x,l) * AHWE(y,l)
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrX)) - LOG10(SUM(PrXY))
enddo

LR = SUM(PrL)

end subroutine Qmerge

! ##############################################################################################

subroutine CheckMerge(SA, SB, kA, kB, LL, focal) 
use Global
implicit none

integer, intent(IN) :: SA, SB, kA, kB, focal
double precision, intent(OUT) :: LL(7)
double precision :: LLtmp(2), ALRtmp(4), LLx(6), LLz(2), LRHS, LLM(6)
integer :: i, j, x, Par(2), SX(2), kx, m, MaybeOpp(2)

LL = 999
ALRtmp = 999
if (kA /= kB .and. focal==1) then
    LL(1) = 777
    return
endif
do i=1, nS(SA, kA)
    if (SibID(i, SA, kA)==GpID(1,SB,kB) .or. SibID(i, SA, kA)==GpID(2,SB,kB)) then
        LL(1) = 777
        exit
    endif
    do j=1, nS(SB, kB)
        if (SibID(j, SB, kB)==GpID(1,SA,kA) .or. SibID(j, SB, kB)==GpID(2,SA,kA)) then
            LL(1) = 777
            exit
        endif
        if (AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB))==999) cycle
        if (AgePriorM(ABS(AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB)))+1,kA) == 0.0) then   
            LL(1) = 777
            exit
        endif
        if (LL(1)==777) exit
    enddo
enddo 
if (LL(1) == 777 .and. focal==1) return

 call Qmerge(SA, SB, kB,  LRHS)
 if (LRHS < -thLR*nS(SA,kA)*nS(SB,kB)) then
    LL(1) = 777
 endif
 if (LL(1) == 777 .and. focal==1) return

 call CalcCLL(SA,kA)
 call CalcCLL(SB,kB)
 call CalcU(-SA,kA, -SB,kB, LL(7))

if (LL(1)/=777 .and. kA==kB) then
    call MergeSibs(SA, SB, kA, LL(1))   ! SB parent of A's
    if (focal==1 .and. (LL(1) > 0 .or. LL(1) - LL(7) < thLR)) return
    call CalcAgeLR(SA, kA, 0, 0, SB, 0, ALRtmp(1))
    if (ALRtmp(1) /=777 .and. LL(1) < 0) then  
        LL(1) = LL(1) + ALRtmp(1)
    else
        LL(1) = 777
    endif
else
    LL(1) = 777
endif
if (focal==1 .and. (LL(1)==777 .or. LL(1) - LL(7) < thLR)) return

 call addFS(0, SA, kA, SB, kB, LL(2))  ! SB FS with an A
 call PairUA(-SB, -SA, kB, kA, LL(3))  ! SB HS with an A
if (LL(2) < 0 .or. LL(3) < 0) then
    call CalcAgeLRCAU(-SB, kB, -SA, kA, ALRtmp(2))
endif
if (ALRtmp(2) /= 777) then
    do x=2,3
        if (LL(x) < 0) then
            LL(x) = LL(x) + ALRtmp(2)
        endif
    enddo
else
    LL(2:3) = 777
endif
 
LLtmp = 999
 call addFS(0, SB, kB, SA, kA, LLtmp(1))   ! SB GP of A's
 call PairUA(-SA, -SB, kA, kB, LLtmp(2))  ! SB GP of A's
 LL(4) = MaxLL(LLtmp)
 call CalcAgeLRCAU(-SA, kA, -SB, kB, ALRtmp(3))
if (ALRtmp(3)/=777) then
    if (LL(4) < 0) then
        LL(4) = LL(4) + ALRtmp(3)
    endif
else
    LL(4) = 777
endif
  
 call ParentHFS(0, SA, kA, SB, kB,3, LL(5))  ! SB FA of A's 
! TODO: PairUA for FS clusters

LLx = 999
 call ParentHFS(0, SA, kA, SB, kB, 1, LLx(1))  !SB HA of A's
 call ParentHFS(0, SA, kA, SB, kB, 2, LLx(2))
if (nAgeClasses > 2) then
    call dummyGP(SA, SB, kA, kB, LLx(3))  ! SB GGP of A's
    call dummyGP(SB, SA, kB, kA, LLx(4))  ! SA GGP of B's
endif
LLz = 999
do x=1,2
    if (GpID(x, SA, kA) > 0) then   ! TODO: more general
        if (Parent(GpID(x, SA, kA), kB)==-SB) then
            LLz(x) = 777
        else
            call PairUA(-SB, GpID(x, SA, kA), kB, 3, LLz(x))
        endif
        if (LLz(x) < 0) then
            LLx(4+x) = LLz(x) + CLL(SA, kA) - Lind(GpID(x, SA, kA))
        endif
    endif
enddo
LL(6) = MaxLL(LLx)  ! most likely 3rd degree relative

LLM = 999
Par = 0
if (kA == kB .and. focal==1 .and. ABS(MaxLL(LL) - LL(1))<thLRrel) then    
! check if they're all FS, in which case merging might also/instead go via 3-k
    SX = (/SA, SB/)
    kx = kA  !=kB
    MaybeOpp = 0
    do m=1,2
        call getFSpar(SX(m), kx, Par(m))
        if (Par(m)>0) cycle
        if (Par(m)==0 .and. ANY(Parent(SibID(1:nS(SX(m),kx),SX(m),kx),3-kx)>0)) cycle
        MaybeOpp(m) = 1
    enddo
    if (MaybeOpp(1)==1 .and. MaybeOpp(2)==1) then   
        if (Par(1)==Par(2) .and. Par(1)/=0) then ! .and. nS(-Par(1),3-kx)==nS(SA,kA)+nS(SB,kB)) then
            MaybeOpp = 0
        else if (Par(1)<0 .and. Par(2)<0) then
            call CalcAgeLR(-Par(1), 3-kA, 0, 0, -Par(2), 0, ALRtmp(4))
            if (ALRtmp(4)==777) MaybeOpp = 0
        endif
    endif
    if (MaybeOpp(1)==1 .and. MaybeOpp(2)==1) then
        call FSMerge(SA,SB,kA, LLM)  ! LLM(1): not, 2:k, 3:3-k, 4:FS
        LLM(2) = MaxLL((/LLM(2), LL(1)-ALRtmp(1)/))  ! merge via k
        LLM(1) = MaxLL((/LLM(1), LL(7)/))  ! do not merge
        if (par(1)<0 .and. par(2)<0) then
            call MergeSibs(-par(1), -par(2), 3-kA, LLM(3))  ! may include indivs not in SA and SB
            call CalcU(Par(1), 3-kA, Par(2), 3-kA, LLM(1))
        endif
        if (MaxLL(LLM)==LLM(4) .and. LLM(4)-LLM(2) > thLR*MAX(nS(SA,kA),nS(SB,kB))) then  
            LL(1) = LLM(4) + ALRtmp(1)  ! FS merge most likely - go ahead.
        else if (LLM(3)<0 .and. LLM(3)-LLM(1) > thLR) then
            if (LLM(3)-LLM(4) > thLRrel) then
                LL(1) = 222 ! more likely that opp. parent need to be merged, or unclear
            else if (Par(1) < 0 .and. Par(2)<0) then
                if (nS(-Par(1),3-kx)+ns(-Par(2),3-kx) > nS(SA,kA)+nS(SB,kB)) then
                    LL(1) = 222   ! LLM(3) and (4) not comparable
                endif
            endif
        endif
    endif
    ! consider merge par, & SA,SB are PO (pairUA). identical LL if complete overlap
    if (LL(1) /= 222 .and.  par(1)<0 .and. par(2)<0) then
        if (nS(SA,kA) + nS(SB,kB) == nS(-Par(1),3-kx)+ns(-Par(2),3-kx)) then  
            if (GpID(kA,SB,kB)==0) then
                GpID(kA,SB,kB) = -SA   ! temp assignment
                call calcCLL(SB, kB)
                call MergeSibs(-par(1), -par(2), 3-kx, LLM(5)) 
                GpID(kA,SB,kB) = 0
                call calcCLL(SB, kB)
                if (LLM(5) - LLM(4) > thLRrel) then 
                    LL(1) = 222
                endif
            endif
            if (GpID(kB,SA,kA)==0) then
                GpID(kB,SA,kA) = -SB   ! temp assignment
                call calcCLL(SA, kA)
                call MergeSibs(-par(1), -par(2), 3-kx, LLM(6)) 
                GpID(kB,SA,kA) = 0
                call calcCLL(SA, kA)
                if (LLM(6) - LLM(4) > thLRrel) then 
                    LL(1) = 222
                endif
            endif
        endif
    endif
endif

end subroutine CheckMerge

! ##############################################################################################

subroutine getFSpar(SA, kA, par)  ! all individuals in SA are FS to eachother
use Global
implicit none

integer, intent(IN) :: SA,  kA
integer, intent(OUT) :: Par
integer :: i, j

Par = 0
do i=1, nS(SA,kA)
    if (Parent(SibID(i,SA,kA), 3-kA)/=0) then
        Par = Parent(SibID(i,SA,kA), 3-kA)
        do j= i+1, nS(SA, kA)
            if (Parent(SibID(j,SA,kA), 3-kA) /= Par .and. &
              Parent(SibID(j,SA,kA), 3-kA)/=0) then
                Par = 0
                return
            endif
        enddo
    endif
enddo

end subroutine getFSpar

! ##############################################################################################

subroutine FSmerge(SA,SB,k, LL)  ! calc LL if SA and SB merged via both pat & mat
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL(4)
integer :: l, x, y, i, u,v, G(2),z
double precision :: PrL(nSnp,4), PrXY(3,3), PrUV(3,3), PrXV(3,3,3,3,4), &
    PrG(3,2), PrX(3), PrTmp(3)

! TODO: currently assumes no grandparents of sibship 3-k, no close inbreeding
! check done in CheckMerge

LL = 999
do i=1,2
    if (GpID(i,SA,k)/=0) then
        if(GpID(i,SA,k)/=GpID(i,SB,k) .and. GpID(i,SB,k)/=0) then
            G(i) = 0  ! shouldn't happen
        else
            G(i) = GpID(i,SA,k)
        endif
    else
        G(i) = GpID(i,SB,k)
    endif
enddo

PrL = 0
do l=1,nSnp 
    do i=1,2
        call ParProb(l, G(i), k, 0,0, PrG(:,i))
    enddo
    do x=1,3
        do z=1,3
            PrTmp(z) = SUM(AKA2P(x,:,z) * PrG(:,1) * PrG(z,2))
        enddo
        PrX(x) = SUM(PrTmp)
    enddo

    do x=1,3  ! P1
        do y=1,3  ! P2
            PrXY(x,y) = 1  ! XPr(2,x,l, sA,k) * AHWE(y,l)
            do i=1,nS(SA,k)
                if (Genos(l, SibID(i,SA,k))/=-9) then
                    PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l, SibID(i,SA,k)), x, y, l)
                endif
            enddo
        enddo
    enddo
    do u=1,3
        do v=1,3
            PrUV(u,v) = 1  ! XPr(2,u,l, sB,k) * AHWE(v,l)
            do i=1,nS(SB,k)
                if (Genos(l, SibID(i,SB,k))/=-9) then
                    PrUV(u,v) = PrUV(u,v) * OKA2P(Genos(l, SibID(i,SB,k)), u, v, l)
                endif
            enddo
        enddo
    enddo

    PrXV = 0
    do x=1,3 
        do y=1,3
            do u=1,3
                do v=1,3
                    PrXV(x,y,u,v,1) = PrXY(x,y) * XPr(2,x,l, sA,k) * AHWE(y,l) * &
                      PrUV(u,v) * XPr(2,u,l, sB,k) * AHWE(v,l)
                    PrXV(x,y,x,v,2) = PrXY(x,y) * PrX(x) * AHWE(y,l) * &
                      PrUV(x,v) * AHWE(v,l)
                enddo
                PrXV(x,y,u,y,3) = PrXY(x,y) * XPr(2,x,l, sA,k) * AHWE(y,l) * &
                  PrUV(u,y) * XPr(2,u,l, sB,k)
            enddo
            PrXV(x,y,x,y,4) = PrXY(x,y) * PrX(x) * AHWE(y,l) * PrUV(x,y)
        enddo
    enddo
    do x=1,4
        PrL(l,x) = LOG10(SUM(PrXV(:,:,:,:,x)))
    enddo
enddo
LL = SUM(PrL,DIM=1)

end subroutine FSmerge

! ##############################################################################################

subroutine MakeFS(A, B)
use Global
implicit none

integer, intent(IN) :: A,B
integer :: x, i, j

i = MIN(A,B)
j = MAX(A,B)
do x=1, nFS(j)   ! u and/or v may already have FS's
    FSID(nFS(i)+x, i) = FSID(x, j)
enddo
nFS(i) = nFS(i) + nFS(j)
nFS(j) = 0

end subroutine MakeFS

! ##############################################################################################

subroutine DoAdd(A, SB, k)
use Global
implicit none

integer, intent(IN) :: A, SB, k
integer :: i, n

if (nS(SB,k) +1 >= maxSibSize) then
    call rexit("Reached maxSibSize ")
endif

do n=1, nS(SB,k)  ! check for FS
    i = SibID(n,SB,k)
    if (nFS(i)==0) cycle
    if (Parent(A, 3-k)/=0 .and. Parent(A, 3-k)==Parent(i, 3-k)) then
        call MakeFS(A, i)
    endif
enddo

SibID(nS(SB,k)+1, SB, k) = A  ! add A to sibship
Parent(A, k) = -SB
nS(SB,k) = nS(SB,k) + 1
 call calcCLL(SB,k)

do n=1,nS(SB,k)   ! update LL of connected sibships
    i = SibID(n,SB,k)
    if (Parent(i,3-k) < 0) then
        call CalcCLL(-Parent(i,3-k), 3-k)
    endif                    
    call CalcLind(i)
enddo
 call calcCLL(SB,k)
 call CalcLind(A)
 
end subroutine DoAdd

! ##############################################################################################

subroutine DoMerge(SA, SB, k)
use Global
implicit none

integer, intent(IN) :: SA, SB, k
integer :: i, j, n, m, x

if (SA/=0) then
    if (nS(SA,k) + nS(SB,k) >= maxSibSize) then
        call rexit("reached maxSibSize")
    endif

    do n=1,nS(SA,k)    ! check for FS
        i = SibID(n,SA,k)
        if (nFS(i)==0) cycle
        do m=1, nS(SB,k)
            j = SibID(m,SB,k)
            if (nFS(j)==0) cycle
            if (Parent(i, 3-k)/=0 .and. Parent(i, 3-k)==Parent(j, 3-k)) then
                call MakeFS(i, j)
            endif
        enddo
    enddo
    do m=1,nS(SB,k)  ! add sibship SB to SA
        SibID(nS(SA,k)+m, SA, k) = SibID(m, SB, k)
        Parent(SibID(m, SB, k), k) = -SA
        do i=1,2
            if (GpID(i, SA, k)==0 .and. GpID(i, SB, k)/=0) then
                GpID(i, SA, k) = GpID(i, SB, k)  ! checked for mismatches earlier
            endif
        enddo
    enddo
    nS(SA,k) = nS(SA,k) + nS(SB,k)
    
    call calcCLL(SA,k)
    do n=1,nS(SA,k) 
        i = SibID(n,SA,k)
        if (Parent(i,3-k) < 0) then
            call CalcCLL(-Parent(i,3-k), 3-k)
        endif                    
        call CalcLind(i)
    enddo
endif
 
do x=SB, nC(k)-1  !remove cluster SB, shift all subsequent ones
    SibID(:, x, k) = SibID(:, x+1, k)
    nS(x, k) = nS(x+1, k)
    GpID(:, x,k) = GpID(:, x+1,k)
    do n=1, nS(x,k)
        Parent(SibID(n,x,k),k) = -x  !Parent(SibID(n,x,k),k)+1 ! shift towards zero.
    enddo
    call calcCLL(x,k)
enddo
SibID(:,nC(k),k) = 0
GpID(:,nC(k),k) = 0
nS(nC(k), k) = 0
nC(k) = nC(k) -1

do x=SB, nC(k)  !fix GPs
    do m=1,2
        do n=1, nC(m)
            if (GpID(k, n, m) == -x)  GpID(k, n, m) = -x+1 
        enddo
    enddo
enddo

end subroutine DoMerge

! ##############################################################################################

subroutine CalcU(A, kA, B, kB, LL)  ! A, SB, k, SA, LL
! TODO: gives incorrect result for A>0, kA==kB, B<0, parent(A,3-kA)=GB
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
double precision, intent(OUT) :: LL
integer :: nA, nB, l,x,y, v, i, j, z,m, f, e, u,n, r, AA(maxSibSize), BB(maxSibSize), &
  PA, PB, DoneA(maxSibSize), catA(maxSibSize), catB(maxSibSize), GA(2), GB(2), par(2),g, Ei
double precision :: PrL(nSnp), PrXY(3,3), PrUZ(3,3, 3,3,3,3, 2), PrPA(3), PrPB(3), &
  PrGA(3,2), PrGB(3,2), PrE(3,2), PrH(3)
 
LL = 999
 
if (A>0) then
    call CalcLind(A)
else if (A<0) then
    call CalcCLL(-A, kA)
endif
if (B>0) then
    call CalcLind(B)
else if (B<0) then
    call CalcCLL(-B, kB)
endif
 
!==================================
if (A==0) then
    if (B==0) then
        LL = 0
    else if (B>0) then
        LL = Lind(B)
    else if (B<0) then
        LL = CLL(-B, kB)
    endif
    return
else if (B==0) then
    if (A>0) then
        LL = Lind(A)
    else if (A<0) then
        LL = CLL(-A,kA)
    endif
    return
endif

if (A>0 .and. B<0) then
    if (Parent(A,kB)==B) then
        LL = CLL(-B,kB)
        return
    endif
else if (B>0 .and. A<0) then
    if (Parent(B,kA)==A) then
        LL = CLL(-A, kA)
        return
    endif
endif

!==================================

PA = 0
PB = 0
if (A>0) then
    call CalcLind(A)
    nA = 1
    AA(1) = A
    PA = Parent(A,kA)
else if (A<0) then
    call CalcCLL(-A, kA)
    nA = nS(-A, kA)
    AA(1:nA) = SibID(1:nA, -A, kA)
    PA = A
endif
if (PA>0) then
    GA = Parent(PA, :)
else if (PA<0) then
    GA = GpiD(:, -PA, kA)
else
    GA = 0
endif

if (B>0) then
    call CalcLind(B)
    nB = 1
    BB(1) = B
    PB = Parent(B,kB)
else if (B<0) then
    call CalcCLL(-B, kB)
    nB = nS(-B, kB)
    BB(1:nB) = SibID(1:nB, -B, kB)
    PB = B
endif
if (PB>0) then
    GB = Parent(PB, :)
else if (PB<0) then
    GB = GpiD(:, -PB, kB)
else
    GB = 0
endif

!================
if (A>0 .and. B>0) then   ! quicker. (?)
    do m=1,2
        if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
           par(m) = Parent(A,m)
        else
            par(m) = 0  ! unknown or unequal
        endif
    enddo
    
     if (Par(1)==0 .and. Par(2)==0) then
        LL = Lind(A) + Lind(B)
        return
    else !if (Par(1)>=0 .and. Par(2)>=0) then
        PrL = 0
        do l=1,nSnp
            if (Genos(l,A)==-9 .and. Genos(l,B)==-9) then
                cycle
            else if (Genos(l,A)==-9) then
                PrL(l) = LindX(l,B)
                cycle
            else if (Genos(l,B)==-9) then
                PrL(l) = LindX(l,A)
                cycle
            endif
 
            if (par(1)/=0 .and. par(2)/=0) then  ! FS
                if (Par(1)>0 .and. Par(2)>0) then
                    call ParProb(l, Par(1), 1, 0,0, PrPA) 
                    call ParProb(l, Par(2), 2, 0,0, PrPB)
                else  
                    do m=1,2
                        if (Par(m) < 0) then
                            if (nS(-Par(m),m)==2) then  ! A + B only members
                                call ParProb(l, Par(m), m, -1,0, PrPA)
                            else
                                call ParProb(l, Par(3-m), 3-m, -1,0, PrH)  ! approx.
                                do x=1,3
                                    if (Genos(l,A)/=-9) then
                                        PrH = PrH * OKA2P(Genos(l,A), x, :, l)
                                    endif
                                    if (Genos(l,B)/=-9) then
                                        PrH = PrH * OKA2P(Genos(l,B), x, :, l)
                                    endif
                                    PrPA(x) = XPr(3, x, l, -Par(m), m)/SUM(PrH)
                                enddo
                                PrPA = PrPA/SUM(PrPA)
                            endif
                            call ParProb(l, Par(3-m), 3-m, 0,0, PrPB)  ! TODO properly: both <0
                        endif
                    enddo
                endif
                do x=1,3
                    do y=1,3
                        PrXY(x,y) = OKA2P(Genos(l,A),x,y,l) * OKA2P(Genos(l,B),x,y,l) * &
                          PrPA(x) * PrPB(y)
                    enddo
                enddo
            else  ! HS
                do m=1,2
                    if (Par(m)==0) cycle
                    call ParProb(l, Parent(A, 3-m), 3-m, A,0, PrPA)
                    call ParProb(l, Parent(B, 3-m), 3-m, B,0, PrPB)
                    call ParProb(l, Par(m), m, A, B, PrH)

                    do x=1,3  ! shared parent
                        do y=1,3  ! parent A
                            if (Parent(A, 3-m) == Parent(B, 3-m) .and. Parent(A, 3-m)/=0) then
                                PrXY(x,y) = OKA2P(Genos(l,A),x,y,l) * PrH(x) * &
                                  PrPA(y) * OKA2P(Genos(l,B),x,y,l)
                            else
                                PrXY(x,y) = OKA2P(Genos(l,A),x,y,l) * PrH(x) * &
                                  PrPA(y) * SUM(OKA2P(Genos(l,B),x,:,l) * PrPB(:))
                            endif
                        enddo
                    enddo
                enddo
            endif
            PrL(l) = LOG10(SUM(PrXY))
        enddo
        LL = SUM(PrL)
        return
    endif
endif

!==================================

if (A>0) then
    if (PA<0) then
        nA = nS(-PA, kA)
        AA(1:nA) = SibID(1:nA, -PA, kA)
    endif
endif
if (B>0 .and. PB<0) then
    nB = nS(-PB, kB)
    BB(1:nB) = SibID(1:nB, -PB, kB)
endif

 catA = 0
 catB = 0
do m = 1, 2
    do i = 1, nA
        do j = 1, nB
            if (AA(i) == BB(j)) cycle
            if (Parent(AA(i), m)<0 .and. Parent(AA(i), m) == Parent(BB(j), m)) then  ! TODO: /=0
                if (m/=kA) then
                    catA(i) = 1
                endif
                if (m/=kB) then
                    catB(j) = 1
                endif
            endif
            if (PA < 0) then
                if (GA(m) < 0) then
!                    if ((A<0 .or. AA(i)==A) .and. GA(m) == Parent(AA(i), m) .and. m/=kA) then
!                        catA(i) = 2
!                    endif
                    if ((B<0 .or. BB(j)==B) .and. GA(m) == Parent(BB(j), m) .and. m/=kB) then
                        catB(j) = 2
                    endif
                endif
            endif
            if (Parent(AA(i), 3-kA) < 0) then
                if (GpID(kB, -Parent(AA(i), 3-kA), 3-kA)==PB .and. PB<0) then
                    catA(i) = 6
                endif
            endif
            if (PB < 0) then
                if (GB(m) < 0) then
                    if ((A<0 .or. AA(i)==A) .and. GB(m) == Parent(AA(i), m) .and. m/=kA) then
                        catA(i) = 3
                    endif
!                    if ((B<0 .or. BB(j)==B) .and. GB(m) == Parent(BB(j), m) .and. m/=kB) then
!                        catB(j) = 3  
!                    endif
                endif
            endif
            if (Parent(BB(j), 3-kB) < 0) then
                if (GpID(kA, -Parent(BB(j),3-kB), 3-kB)==PA .and. PA<0) then
                    catB(j) = 6
                endif  ! TODO: not considered yet GpID(kA, -Parent(BB(j),3-kB), 3-kB) == Parent(AA(i), 3-kA)
            endif
        enddo
    enddo
enddo
if (PB<0 .and. PB==GA(kB) .and. PA<=0) then
    catA(nA+1) = 4
endif  
if (PA<0 .and. PA==GB(kA) .and. PB<=0) then
    catB(nB+1) = 4 
endif
do n = 1, 2
    if (GA(n) < 0 .and. PA<0 .and. PB<0) then 
        if (GA(n) == GB(n)) then
            catA(nA+1) = 5  
            catB(nB+1) = 5
        endif
    endif
enddo
if (kA==kB .and. PA==PB .and. PA/=0) then
    catA(nA+1) = 7 
    catB(nB+1) = 7
endif

if (ALL(catA==0) .and. ALL(catB==0)) then
    if (A>0 .and. B>0) then
        LL = Lind(A) + Lind(B)
    else if (A<0 .and.B>0) then
        LL = CLL(-A,kA) + Lind(B)
    else if (A>0 .and. B<0) then
        LL = Lind(A) + CLL(-B,kB)
    else if (A<0 .and. B<0) then
        LL = CLL(-A,kA) + CLL(-B,kB)
    endif
    
    return
endif
  
LL = 999
PrL = 0
do l=1, nSnp
    PrUZ = 0
    ! == grandparents ==
    if (PA>0) then   !    .and. catA(nA+1)/=4 not: PA member of SB
        call ParProb(l, PA, kA, 0,0, PrPA)
    else
        PrPA = 1
    endif
    if (PB>0) then  ! .and. catB(nB+1)/=5
        call ParProb(l, PB, kB, 0,0, PrPB)
    else
        PrPB = 1
    endif
    do m=1, 2
        if ((ANY(catA==2) .and. m/=kA) .or. (ANY(catB==2) .and. m/=kB)) then  ! catA(nA+1)==5
            call ParProb(l, GA(m), m, -1,0, PrGA(:, m))
        else if (catA(nA+1)==4 .and. m/=kB .and. ANY(Parent(BB(1:nB), 3-kB)==GA(m) .and. GA(m)<0)) then
            call ParProb(l, GA(m), m, -1,0, PrGA(:, m))
            do e=1,3
                do j=1, nS(-GA(m), m)
                    Ei = SibID(j, -GA(m), m)
                    if (nFS(Ei) == 0) cycle
                    if (Parent(Ei, 3-m) == GA(3-m) .and. GA(3-m)/=0) cycle    ! should be redundant
                    if (Parent(Ei, kB) == PB .and. PB/=0) cycle  
                    if (Parent(Ei, kA) == PA .and. PA/=0) cycle  
                    call ParProb(l, Parent(Ei, 3-m), 3-m, Ei,-1, PrH) 
                    do i=1, nFS(Ei)
                        if (Genos(l,FSID(i, Ei))==-9) cycle
                        if (FSID(i,Ei) == A) cycle
                        if (FSID(i,Ei) == B) cycle
                        PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e, l)
                    enddo
                    PrGA(e, m) = PrGA(e, m) * SUM(PrH)
                enddo
            enddo
            PrGA(:, m) = PrGA(:, m) / SUM(PrGA(:,m))
        else
            call ParProb(l, GA(m), m, 0,0, PrGA(:, m))
        endif
        if ((ANY(catA==3) .and. m/=kA) .or. (ANY(catB==3) .and. m/=kB)) then  ! PB>0 .or. 
            call ParProb(l, GB(m), m, -1,0, PrGB(:, m))
        else if (catB(nB+1)==4 .and. m/=kA .and. ANY(Parent(AA(1:nA), 3-kA)==GB(m) .and. GB(m)<0)) then  ! optional for catA(nA+1)==4
            call ParProb(l, GB(m), m, -1,0, PrGB(:, m))
            do e=1,3
                do j=1, nS(-GB(m), m)
                    Ei = SibID(j, -GB(m), m)
                    if (nFS(Ei) == 0) cycle
                    if (Parent(Ei, 3-m) == GB(3-m) .and. GB(3-m)/=0) cycle    ! should be redundant
                    if (Parent(Ei, kB) == PB .and. PB/=0) cycle  
                    if (Parent(Ei, kA) == PA .and. PA/=0) cycle  
                    call ParProb(l, Parent(Ei, 3-m), 3-m, Ei,-1, PrH) 
                    do i=1, nFS(Ei)
                        if (Genos(l,FSID(i, Ei))==-9) cycle
                        if (FSID(i,Ei) == A) cycle
                        if (FSID(i,Ei) == B) cycle
                        PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e, l)
                    enddo
                    PrGB(e, m) = PrGB(e, m) * SUM(PrH)
                enddo
            enddo
            PrGB(:, m) = PrGB(:, m) / SUM(PrGB(:,m))
        else
            call ParProb(l, GB(m), m, 0,0, PrGB(:, m))
        endif    
    enddo

    do x=1,3 
        do y=1,3
            do u=1,3  ! GP A, kB
                do z=1,3  ! GP A, 3-kB
                    do v=1,3  ! GP B, kB
                        if (kA==kB .and. PA==PB .and. PA/=0) then
                            PrUZ(x,x,u,z,u,z,1) = PrPA(x) * PrPB(x) * &
                              AKA2P(x,u,z) * PrGA(u,kB) * PrGA(z,3-kB)
                        else if (catA(nA+1)==5 .or. catB(nB+1)==5) then
                            if (GA(1)/=0 .and. GA(1) == GB(1) .and. GA(2)/=0 .and. GA(2)==GB(2)) then
                                PrUZ(x,y,u,z,u,z,1) = PrPA(x) * PrPB(y) * &
                                  AKA2P(x,u,z) * AKA2P(y,u,z) * PrGA(u,kB) * PrGA(z,3-kB)
                            else
                                do m=1,2
                                    if (GA(m)/=0 .and. GA(m) == GB(m)) then
                                        if (m==kB) then
                                            PrUZ(x,y,u,z,u,:,1) = PrPA(x) * PrPB(y)* AKA2P(x,u,z) * &
                                              AKA2P(y,u,:) * PrGA(u,m) * PrGA(z,3-m) * PrGB(:, 3-m)   
                                        else
                                            PrUZ(x,y,u,z,v,z,1) = PrPA(x) * PrPB(y) * AKA2P(x,u,z) * &
                                              AKA2P(y,v,z) * PrGA(u,3-m) * PrGA(z,m) * PrGB(v, 3-m)
                                        endif
                                    endif
                                enddo
                            endif
                        else if (catA(nA+1)==4) then   ! PB = GA(kB)
                            if (kA==kB .and. GA(3-kB)==GB(3-kB) .and. GA(3-kB)/=0) then  ! inbreeding loop
                                PrUZ(x,y,y,z,v,z,1) = PrPA(x) * PrPB(y) * AKA2P(x,y,z) * &
                                  AKA2P(y,v,z) * PrGA(z,3-kB) * PrGB(v,kB)
                            else
                                PrUZ(x,y,y,z,v,:,1) = PrPA(x) * PrPB(y) * AKA2P(x,y,z) * &
                                AKA2P(y,v,:) * PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:, 3-kB)   
                            endif
                        else if (catB(nB+1)==4) then  ! PA = GB(kA)
                            if (kA==kB) then
                                if (GA(3-kA)==GB(3-kA) .and. GA(3-kA)/=0) then  ! inbreeding loop
                                    PrUZ(x,y,u,z,x,z,1) = PrPA(x) * PrPB(y) * AKA2P(x,u,z) * &
                                      AKA2P(y,x,z) * PrGA(z,3-kB) * PrGA(u,kB)
                                else
                                    PrUZ(x,y,u,z,x,:,1) = PrPA(x) * PrPB(y) * AKA2P(x,u,z) * &
                                      AKA2P(y,x,:) * PrGA(z,3-kB) * PrGA(u,kB) * PrGB(:, 3-kB)   
                                endif
                            else if (kA/=kB) then
                                if (GA(3-kA)==GB(3-kA) .and. GA(3-kA)/=0) then  ! inbreeding loop
                                    PrUZ(x,y,u,z,u,x,1) = PrPA(x) * PrPB(y) * AKA2P(x,u,z) * &
                                      AKA2P(y,u,x) * PrGA(z,3-kB) * PrGA(u,kB)
                                else
                                    PrUZ(x,y,u,z,v,x,1) = PrPA(x) * PrPB(y) * AKA2P(x,u,z) * &
                                      AKA2P(y,v,x) * PrGA(z,3-kB) * PrGA(u,kB) * PrGB(v,kB)
                                endif
                            endif
                        else
                            PrUZ(x,y,u,z,v,:,1) = PrPA(x) * PrPB(y) * &  
                              AKA2P(x,u,z) * AKA2P(y,v,:) * &
                              PrGA(u,kB) * PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:, 3-kB)
                        endif
                    enddo
                enddo
            enddo
            PrXY(x,y) = SUM(PrUZ(x,y, :,:,:,:, 1))
            PrUZ(x,y, :,:,:,:, 2) = PrUZ(x,y, :,:,:,:, 1)
        enddo
    enddo
    
    ! == siblings ==

    do x=1,3  ! SA
        doneA = 0
        do y=1,3  ! SB
            do j=1, nB
                if (nB>1 .and. nFS(BB(j))==0) cycle
                if ((catB(j)==1 .and. kA/=kB) .or. catB(j)==2 .or. catB(j)==3) then
                    PrE = 1
                else if (catB(j)==6) then
                    call ParProb(l, GpID(3-kA, -Parent(BB(j), 3-kB), 3-kB), 3-kA, 0,0, PrH)
                    do e=1,3
                        PrE(e,1) = SUM(AKA2P(e,x,:) * PrH) 
                    enddo
                else if (B<0 .or. nB==1 .or. (nB>1 .and. ANY(FSID(1:nFS(BB(j)), BB(j))==B))) then
                    call ParProb(l, Parent(BB(j),3-kB), 3-kB, -1,0, PrE(:,1))
                else
                    call ParProb(l, Parent(BB(j),3-kB), 3-kB, BB(j),-1, PrE(:,1))
                endif
                PrE(:,2) = PrE(:,1)  ! 1=all; 2=non-sibs only
                
                do e=1,3
                    do f=1, MAX(nFS(BB(j)),1)  ! includes some AA if kA/=kB & cat=1 & A<0
                        if (B<0 .and. nFS(BB(j))==0) cycle
                        if (Genos(l,FSID(f, BB(j)))==-9) cycle
                        if (B<0 .or. (nB==1 .and. FSID(1, BB(j))==B) .or. (nB>1 .and. FSID(f, BB(j))==B)) then 
                            PrE(e,1) = PrE(e,1) * OKA2P(Genos(l,FSID(f, BB(j))), y, e, l)
                        else if (nB>1 .and. .not. (kA/=kB .and. Parent(BB(j),3-kB)==PA .and. PA/=0)) then    ! TODO double check
                           PrE(e,:) = PrE(e,:) * OKA2P(Genos(l,FSID(f, BB(j))), y, e, l)
                        endif
                    enddo

                    if (catB(j)==1 .and. kA==kB) then  !  Parent(BB(j), 3-kB)/=0, shared opposite parent
                        do i=1,nA
                            if (nA>1 .and. nFS(AA(i))==0) cycle
                            if (Parent(AA(i), 3-kB) /= Parent(BB(j), 3-kB)) cycle
                            do f=1, MAX(nFS(AA(i)),1)
                                if (A<0 .and. nFS(AA(i))==0) cycle
                                if (Genos(l,FSID(f, AA(i)))==-9) cycle
                                if (A<0 .or. (nA==1 .and. FSID(1,AA(i))==A) .or. (nA>1 .and. FSID(f,AA(i))==A)) then  
                                    PrE(e,1) = PrE(e,1) * OKA2P(Genos(l,FSID(f,AA(i))), x, e, l)
                                    DoneA(i) = 1
                                else if (nA>1) then
                                    PrE(e,:) = PrE(e,:) * OKA2P(Genos(l,FSID(f,AA(i))), x, e, l)
                                    DoneA(i) = 10
                                endif
                            enddo
                        enddo
                    endif
                
                    if (Parent(BB(j),3-kB) < 0 .and. (B<0 .or. nB==1 .or. &
                      (nB>1 .and. ANY(FSID(1:nFS(BB(j)), BB(j))==B)))) then 
                        do g=1, nS(-Parent(BB(j),3-kB), 3-kB)
                            Ei = SibID(g, -Parent(BB(j),3-kB), 3-kB)
                            if (nFS(Ei) == 0) cycle
                            if (catB(j)==1 .and. kA/=kB) cycle  ! Parent(BB(j),3-kB) == PA
                            if (nB>1 .and. Parent(Ei, kB) == PB .and. PB/=0) cycle 
                            if (Parent(Ei, kA) == PA .and. PA/=0) cycle 
                            if (nB==1 .and. Parent(Ei, kB)==PB .and. PB/=0) then
                                do i=1, nFS(Ei)                  
                                    if (Genos(l,FSID(i, Ei))==-9) cycle
                                    if (FSID(i,Ei) == A) cycle
                                    if (FSID(i,Ei) == B) cycle
                                    PrE(e,:) = PrE(e,:) * OKA2P(Genos(l,FSID(i,Ei)), y, e, l)
                                enddo
                            else
                                call ParProb(l, Parent(Ei, kB), kB, Ei,-1, PrH) 
                                do i=1, nFS(Ei)
                                    if (Genos(l,FSID(i, Ei))==-9) cycle
                                    if (FSID(i,Ei) == A) cycle
                                    if (FSID(i,Ei) == B) cycle
                                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e, l)
                                enddo
                                PrE(e,:) = PrE(e,:) * SUM(PrH)
                            endif
                        enddo
                    endif
                enddo
                
                do r=1,2
                    if (catB(j)==1 .and. kA/=kB) then  ! Parent(BB(j), 3-kB)==PA
                        PrUZ(x,y,:,:,:,:,r) = PrUZ(x,y,:,:,:,:,r) * PrE(x,r)
                    else if (CatB(j)==2) then
                        do z=1,3
                            PrUZ(x,y,:,z,:,:,r) = PrUZ(x,y,:,z,:,:,r) * PrE(z,r)   
                        enddo
                    else if (CatB(j)==3) then
                        do z=1,3
                            PrUZ(x,y,:,:,:,z,r) = PrUZ(x,y,:,:,:,z,r) * PrE(z,r)   
                        enddo   
                    else 
                        PrUZ(x,y,:,:,:,:,r) = PrUZ(x,y,:,:,:,:,r) * SUM(PrE(:,r))
                    endif
                enddo
            enddo
 !       enddo  ! y
        
            do i=1, nA
                if (DoneA(i)==1) cycle
                if (nA>1 .and. nFS(AA(i))==0) cycle
                if (ANY(BB(1:nB) == AA(i)) .and. B<0) cycle  !   kA/=kB .and. Parent(AA(i), 3-kA)==PB .and. PB/=0
                if (catA(i)>0 .and. catA(i)<4) then
                    PrE = 1
                else if (catA(i)==6) then ! Parent(AA(i), 3-kA) <0
                    call ParProb(l, GpID(3-kB, -Parent(AA(i), 3-kA), 3-kA), 3-kB, 0,0, PrH)
                    do e=1,3
                        PrE(e,1) = SUM(AKA2P(e,y,:) * PrH) 
                    enddo
                 else if (A<0 .or. nA==1 .or. (nA>1 .and. ANY(FSID(1:nFS(AA(i)),AA(i))==A))) then  
                    call ParProb(l, Parent(AA(i), 3-kA), 3-kA, -1,0, PrE(:,1))   
                else
                    call ParProb(l, Parent(AA(i), 3-kA), 3-kA, AA(i),-1, PrE(:,1)) 
                endif
                PrE(:,2) = PrE(:,1)  ! 1=all; 2=non-sibs only
                
                do e=1,3
                    do f=1, MAX(nFS(AA(i)),1)
                        if (Genos(l,FSID(f, AA(i)))==-9) cycle
                        if (A<0 .and. nFS(AA(i))==0) cycle
                        if (A<0 .or. nA==1 .or. (nA>1 .and. FSID(f, AA(i))==A)) then                  
                            PrE(e,1) = PrE(e,1) * OKA2P(Genos(l,FSID(f,AA(i))), x, e, l)
                        else if (nA>1 .and. .not. (kA/=kB .and. Parent(AA(i),3-kA)==PB .and. PB/=0)) then  ! done under B
                            PrE(e,:) = PrE(e,:) * OKA2P(Genos(l,FSID(f,AA(i))), x, e, l) 
                        endif
                    enddo
                                
                    if (Parent(AA(i), 3-kA) < 0 .and. (A<0 .or. nA==1 .or. &
                     (nA>1 .and. ANY(FSID(1:nFS(AA(i)),AA(i))==A)))) then   
                        do g=1, nS(-Parent(AA(i), 3-kA), 3-kA)
                            Ei = SibID(g, -Parent(AA(i), 3-kA), 3-kA)
                            if (nFS(Ei) == 0) cycle
                            if (nA>1 .and. Parent(Ei, kA) == PA .and. PA/=0) cycle  ! .and. A<0
                            if (Parent(Ei, kB) == PB .and. PB/=0) cycle  !  .and. nFS(B)>1
                            if (nA==1 .and. Parent(Ei, kA)==PA .and. PA/=0) then
                                do j=1, nFS(Ei)                  
                                    if (Genos(l,FSID(j, Ei))==-9) cycle
                                    if (FSID(j,Ei) == A) cycle
                                    if (FSID(j,Ei) == B) cycle
                                    PrE(e,:) = PrE(e,:) * OKA2P(Genos(l,FSID(j,Ei)), x, e, l)
                                enddo
                            else
                                call ParProb(l, Parent(Ei, kA), kA, Ei,-1, PrH)  ! Assume is not one of GPs
                                do j=1, nFS(Ei)                  
                                    if (Genos(l,FSID(j, Ei))==-9) cycle
                                    if (FSID(j,Ei) == A) cycle
                                    if (FSID(j,Ei) == B) cycle
                                    PrH = PrH * OKA2P(Genos(l,FSID(j,Ei)), :, e, l)
                                enddo
                                PrE(e,:) = PrE(e,:) * SUM(PrH)
                            endif
                        enddo
                    endif
                enddo
                
                do r=1,2
                    if (catA(i)==1 .and. kA/=kB) then
                        PrUZ(x,y,:,:,:,:,r) = PrUZ(x,y,:,:,:,:,r) * PrE(y,r)
                    else if (catA(i)==2) then
                        do z=1,3
                            if (kA==kB) then
                                PrUZ(x,y,:,z,:,:,r) = PrUZ(x,y,:,z,:,:,r) * PrE(z,r)
                            else
                                PrUZ(x,y,z,:,:,:,r) = PrUZ(x,y,z,:,:,:,r) * PrE(z,r)         
                            endif
                        enddo
                    else if (catA(i)==3) then  !  .or. (catA(nA+1)==4 .and. Parent(AA(i), 3-kA) == GB(3-kA) .and. GB(3-kA)/=0)
                        do z=1,3
                            if (kA==kB) then
                                PrUZ(x,y,:,:,:,z,r) = PrUZ(x,y,:,:,:,z,r) * PrE(z,r)
                            else
                                PrUZ(x,y,:,:,z,:,r) = PrUZ(x,y,:,:,z,:,r) * PrE(z,r)
                            endif
                        enddo
                    else
                        PrUZ(x,y,:,:,:,:,r) = PrUZ(x,y,:,:,:,:,r) * SUM(PrE(:,r))    
                    endif
                enddo
            enddo  ! i
        enddo  ! x
    enddo  ! y
    PrL(l) = LOG10(SUM(PrUZ(:,:, :,:,:,:, 1)) / SUM(PrUZ(:,:, :,:,:,:, 2)))
enddo

LL = SUM(PrL)


end subroutine CalcU

! ##############################################################################################

subroutine AddSib(A, SB, k, LL)  ! TODO?: don't assume 3-k's remain fixed
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, Bj, m, f, AncB(2,16), Inbr
double precision :: PrL(nSnp), PrX(3), PrY(3)

LL = 999
if (Parent(A,k)==-SB) then
    LL = 888
else if (Parent(A,k)/=0) then
    LL = 777
endif
do f=1, nS(SB,k)
    Bj = SibID(f, SB, k)
    if (Parent(A, 3-k) < 0) then  ! TODO: /=0
        if (Parent(Bj, 3-k) == Parent(A, 3-k)) then
            LL = 777  ! use addFS() instead
        endif
    endif
    if (AgeDiff(A,Bj)==999) cycle
    if (AgePriorM(ABS(AgeDiff(A, Bj))+1,k)==0.0) then   
        LL=777
    endif 
enddo
if (LL/=999) return

 call GetAncest(-SB, k, AncB)
if (ANY(AncB == A)) then  ! A>0
    LL = 777
else
    do x=3, 16
        do m=1,2
            if (AncB(m, x) > 0) then
                if (AgeDiff(A, AncB(m, x)) <=0) then  ! A older than putative ancestor
                    LL = 777
                endif
            endif
        enddo
    enddo
endif

Inbr = 0
if (Parent(A,3-k) < 0) then
    do m = 1,2
        if (GpID(m, -Parent(A,3-k), 3-k) /=0) then
            if (GpID(m, -Parent(A,3-k), 3-k) /= GpID(m, SB, 3-k) .and. &
            GpID(m, SB, 3-k)/=0) then
                LL = 777
           endif
        endif
    enddo
    if (Parent(A,3-k) == GpID(3-k, SB, k)) then
        Inbr = 1  ! inbreeding loop created
    endif
endif
if (LL == 777) return

if (Inbr == 0) then
    PrL = 0
    do l=1,nSnp
        call ParProb(l, Parent(A,3-k), 3-k, A,0, PrY)
        do x=1,3
            PrX(x) = XPr(3,x,l, SB,k)
            if (Genos(l,A) /= -9) then
                PrX(x) = PrX(x) * SUM(OKA2P(Genos(l,A), x, :, l) * PrY)
            endif
        enddo
        PrL(l) = LOG10(SUM(PrX))   
    enddo
    LL = SUM(PrL)
else
    call DoAdd(A, SB, k)
    call CalcU(A, k, -SB, k, LL)
    call RemoveSib(A, SB, k)
endif

end subroutine AddSib

! ##############################################################################################

subroutine MergeSibs(SA, SB, k, LL)  
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL
integer :: l,x,y, r,v, Bj, j, Ai, i, DoneA(nS(SA,k)), G(2), m, Ei, w, z, &
  doneE(nInd), AncA(2,16), AncB(2,16)
double precision :: PrL(nSnp,2), PrX(3), PrY(3), PrG(3,2), PrXV(3,3,3), &
LLUX, PrXW(3,3), PrM(3,3,3), PrZ(3)
logical :: SharedPar

LL = 999
G = 0  
do m=1,2  
    if (GpID(m,SA,k) /= 0) then
        if (GpID(m,SB,k) /= 0 .and. GpID(m,SA,k)/=GpID(m,SB,k)) then
            LL = 777  ! incompatible grandparents
        else
            G(m) = GpID(m,SA,k)  ! including if GP(B) is dummy
        endif
    else if (GpID(m,SA,k) == 0) then
        G(m) = GpID(m,SB,k)
    endif
enddo
if (GpID(k, SA,k)==-SB .or. GpID(k, SB, k)==-SA) then
    LL = 777
endif
if (LL==777) return
 call GetAncest(-SA, k, AncA)
 call GetAncest(-SB, k, AncB)
if (ANY(AncA(k, 3:16) == -SB)) LL = 777
if (ANY(AncB(k, 3:16) == -SA)) LL = 777
if (LL==777) return 

! are the sibships currently connected?
sharedPar = .FALSE.
do r = 1, nS(SB, k)
    Bj = SibID(r, SB, k)  
    if (NFS(Bj) == 0) cycle
    if (Parent(Bj, 3-k) == 0) cycle
    do v = 1, nS(SA,k)  ! condition on parent(A,3-k)
        Ai = SibID(v, SA, k)
        if (nFS(Ai) == 0) cycle
        if (Parent(Ai, 3-k)==Parent(Bj,3-k) .or. &
          (Parent(Ai, 3-k) <0 .and. Parent(Bj, 3-k)<0)) then
            sharedPar = .TRUE.
            exit    
        endif
    enddo
enddo

PrL = 0
if (.NOT. sharedPar) then
    do l=1,nSnp
        do m=1,2
            call ParProb(l, G(m), m, 0,0, PrG(:,m))
        enddo  
        PrXV = 0
        do r=1,3  ! G1
            do v=1,3  ! G2
                do x=1,3
                    PrXV(x,r,v) = XPr(1,x,l, SA,k) * XPr(1,x,l, SB,k) * &
                        AKA2P(x,r,v) * PrG(r,1) * PrG(v,2) 
                enddo
            enddo
        enddo
        PrL(l,1) = LOG10(SUM(PrXV))
    enddo
    LL = SUM(PrL(:,1))

else if (sharedPar) then
    do l=1,nSnp
        if (G(1)==0 .and. G(2)==0) then
            PrX = AHWE(:,l)
        else
            do m=1,2
                call ParProb(l, G(m), m, -1,0, PrG(:,m))
            enddo
            do x=1,3
                do y=1,3
                    PrY(y) = SUM(AKA2P(x, y, :) * PrG(y,1) * PrG(:,2))  ! grandparents, when merged
                enddo
                PrX(x) = SUM(PrY)
            enddo
        endif
        do x=1,3
            do w=1,3
                PrXW(x,w) = XPr(2,x,l,SB,k) * XPr(2,w,l,SA,k)  ! gps, when separate
            enddo
        enddo
        
        DoneA = 0
        DoneE = 0
        do r=1, nS(SB,k)
            Bj = SibID(r, SB, k) 
            if (NFS(Bj) == 0) cycle  ! moved to its FS
            call ParProb(l, Parent(Bj,3-k), 3-k, -1,0, PrY)
            do y=1,3   
                PrM(:,:,y) = PrY(y) 
                do x=1,3 
                    do j=1, nFS(Bj)  ! default: nFS = 1
                        if (Genos(l,FSID(j, Bj))/=-9) then
                            PrM(x,:,y) =  PrM(x,:,y) * OKA2P(Genos(l,FSID(j,Bj)), x, y, l)
                        endif
                    enddo
                enddo

                if (Parent(Bj,3-k) < 0) then
                    do v = 1, nS(-Parent(Bj, 3-k), 3-k)
                        Ei = SibID(v, -Parent(Bj, 3-k), 3-k)
                        if (Parent(Ei, k) == -SA .or. Parent(Ei, k) == -SB) cycle
                        if (nFS(Ei) == 0) cycle
                        call ParProb(l, Parent(Ei, k), k, Ei,-1, PrZ)                      
                        do z=1,3
                            do i=1, nFS(Ei)
                                if (Genos(l,FSID(i, Ei))/=-9) then
                                    PrZ(z) = PrZ(z) * OKA2P(Genos(l,FSID(i,Ei)), z, y, l)
                                endif
                            enddo
                        enddo
                        PrM(:,:,y) = PrM(:,:,y) * SUM(PrZ)
                        DoneE(SibID(v, -Parent(Bj, 3-k), 3-k)) = 1
                    enddo
                endif
                
                if (Parent(Bj,3-k)/=0) then
                    do v = 1, nS(SA,k)
                        Ai = SibID(v, SA, k)
                        if (nFS(Ai) == 0) cycle
                        if (Parent(Ai, 3-k)/=Parent(Bj,3-k)) cycle
                        do w=1,3
                            do i=1, nFS(Ai)
                                if (Genos(l,FSID(i, Ai))/=-9) then
                                    PrM(:,w,y) =  PrM(:,w,y) * OKA2P(Genos(l,FSID(i,Ai)), w, y, l)
                                endif
                            enddo
                        enddo
                        doneA(v) = 1
                    enddo
                endif
            enddo ! y
            do x=1,3
                PrX(x) = PrX(x) * SUM(PrM(x,x,:))
                do w=1,3
                    PrXW(x,w) = PrXW(x,w) * SUM(PrM(x,w,:))
                enddo
            enddo
        enddo  ! r=1, nS(SB,k)
            
        do v = 1, nS(SA,k)  
            if (doneA(v)==1) cycle
            Ai = SibID(v, SA, k)
            if (NFS(Ai) == 0) cycle
            call ParProb(l, Parent(Ai,3-k), 3-k, -1,0, PrY)
            do y=1,3
                PrM(:,:,y) = PrY(y)
                do w=1,3
                    do i=1, nFS(Ai)
                        if (Genos(l,FSID(i, Ai))/=-9) then
                            PrM(:,w,y) =  PrM(:,w,y) * OKA2P(Genos(l,FSID(i,Ai)), w, y, l)
                        endif
                    enddo
                enddo
                
                if (Parent(Ai,3-k) < 0) then
                    do r = 1, nS(-Parent(Ai, 3-k), 3-k)
                        Ei = SibID(r, -Parent(Ai, 3-k), 3-k)
                        if (doneE(Ei)==1) cycle
                        if (Parent(Ei, k) == -SA .or. Parent(Ei, k) == -SB) cycle
                        if (nFS(Ei) == 0) cycle
                        call ParProb(l, Parent(Ei, k), k, Ei,-1, PrZ)
                        do z=1,3
                            do i=1, nFS(Ei)
                                if (Genos(l,FSID(i, Ei))/=-9) then
                                    PrZ(z) = PrZ(z) * OKA2P(Genos(l,FSID(i,Ei)), z, y, l)
                                endif
                            enddo
                        enddo
                        PrM(:,:,y) =  PrM(:,:,y) * SUM(PrZ)
                    enddo
                endif
            enddo
            do x=1,3
                PrX(x) = PrX(x) * SUM(PrM(x,x,:))
                do w=1,3
                    PrXW(x,w) = PrXW(x,w) * SUM(PrM(x,w,:))
                enddo
            enddo
        enddo

        PrL(l,1) = LOG10(SUM(PrX))  ! merged
        PrL(l,2) = LOG10(SUM(PrXW))  ! cond. unrelated
    enddo
    
    call CalcU(-SA,k, -SB, k, LLUX)   
    LL = SUM(PrL(:,1)) - SUM(PrL(:,2)) + LLUX
endif

end subroutine MergeSibs

! ##############################################################################################

subroutine AddFS(A, SB, kB, SA, kA, LL)  ! A/SA FS with any B?
use Global
implicit none

integer, intent(IN) :: A, SB, kB, SA, kA
double precision, intent(OUT) :: LL
integer :: l, x, y, Par(nS(SB,kB)), i, Bj, Ei, f, g, MaybeFS(nS(SB,kB)), z, &
 PA, AncA(2,16), AncB(2,16), AncPF(2,16), m, h
double precision :: PrL(nSnp, nS(SB,kB),2), PrY(3,2), PrX(3,2), PrZ(3), &
  dLL(nS(SB,kB)), LLtmp(2), LLUX

PrL = 0
LL = 999

Par = 0  ! shared parent 3-kB  (cand. parent(kB) == SB)
MaybeFS = 1

 call GetAncest(-SB, kB , AncB)
PA = 0
if (A /= 0) then
    PA = Parent(A, 3-kB)
    call GetAncest(A, kA, AncA)
else if (SA /= 0) then   ! TODO: does it matter if kA=kB?
    PA = GpID(3-kB, SA, kA)
     call GetAncest(-SA, kA, AncA)
endif

if (A/=0) then
    if (Parent(A,kB)/=0 .and. Parent(A,kB)/=-SB) then
        LL = 777
    else if (ANY(AncB == A)) then 
        LL = 777
    else
        do x=2, 16
            do m=1,2
                if (AncB(m, x) > 0) then
                    if (AgeDiff(A, AncB(m, x)) <=0) then  ! A older than putative ancestor
                        LL = 777
                    endif
                endif
            enddo
        enddo
    endif
else if (SA/=0) then
if (GpID(kB, SA, kA)/=0 .and. GpID(kB, SA, kA)/=-SB) then
        LL = 777
    else if (ANY(AncB(kA, 2:16) == -SA)) then
        LL = 777
    else
        do x=2, 16
            do m=1,2
                if (AncB(m, x) > 0) then
                    if (ANY(AgeDiff(SibID(1:nS(SA,kA),SA,kA), AncB(m, x)) <=0)) then  ! A older than putative ancestor
                        LL = 777
                    endif
                endif
            enddo
        enddo
    endif
endif
if (LL /= 999) return

if (ANY(AncA(kB, 3:16) == -SB)) then
    LL = 444   ! TODO: check
    return
endif

do f=1, nS(SB,kB)
    if (NFS(SibID(f, SB, kB))==0) then
        MaybeFS(f) = -1
        cycle
    endif   
    do i=1,nFS(SibID(f, SB, kB))
        Bj = FSID(i, SibID(f, SB, kB))
        if (A == Bj) then
            LL = 888
        else if (A >0) then
            if (Parent(A,3-kB) == Bj) then  
               MaybeFS(f) = 0          ! possible, but unlikely    
            else if (Parent(Bj, 3-kB) == A) then
                MaybeFS(f) = 0
            else if (AgeDiff(A,Bj)/=999) then
                if (AgePriorM(ABS(AgeDiff(A, Bj))+1, kB)==0.0) then   
                    LL=777
                else if (AgePriorM(ABS(AgeDiff(A, Bj))+1,3-kB)==0.0) then   
                    MaybeFS(f) = 0
                endif
            endif
        else if (SA/=0 .and. kA/=kB) then
            if (Parent(Bj, 3-kB) == -SA) then
                MaybeFS(f) = 0  ! cannot be FS with own parent
            endif
        endif
        if (Bj == PA .or. (A/=0 .and. A == Parent(Bj, 3-kB))) then
            MaybeFS(f) = 0
            cycle
        endif
        if (PA>0) then
            if (Parent(PA,1)==Bj .or. Parent(PA,2)==Bj) then
                MaybeFS(f) = 0
                cycle
            endif
        endif
        
        Par(f) = Parent(Bj, 3-kB)
        if (PA/=0 .and. PA/=Par(f) .and. Par(f)/=0) then
            MaybeFS(f) = 0
        else if (Par(f)==0) then
            Par(f) = PA
        endif
    enddo
    if (Par(f)<0 .and. SA/=0 .and. kA==kB) then
        do i=1, nS(SA,kA)
            if (Parent(SibID(i,SA,kA), 3-kA) == Par(f)) then
                LL = 444   ! TODO: implement (P-O inbreeding)
            endif
        enddo
    endif
enddo
if (LL /= 999) return
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
    LL = 777
    return
endif

do f=1, nS(SB,kB)
    if (nFS(SibID(f, SB, kB))==0) cycle
    if (MaybeFS(f)/=1 .or. Par(f)==0 .or. Par(f)==PA) cycle
    call getAncest(Par(f), 3-kB, AncPF)
    if (Par(f)>0) then
        if (A/=0 .and. ANY(AncPF(:,2:16)==A)) MaybeFS(f) = 0
        if (SA/=0 .and. ANY(AncPF(kA,2:16)==-SA)) MaybeFS(f) = 0
    else if (Par(f) < 0) then
        if (A/=0 .and. ANY(AncPF(:,3:16)==A)) MaybeFS(f) = 0
        if (SA/=0 .and. ANY(AncPF(kA,3:16)==-SA)) MaybeFS(f) = 0
        if (SA/=0 .and. GpID(kA,SB,kB)==0 .and. GpID(kA,-Par(f),3-kB)==0) then
            call PairUA(SibID(f,SB,kB),-SA,kB, kA, LLtmp(1))
            if (LLtmp(1) < 0) then ! can't tel if FS or double GP
                MaybeFS(f) = 0
            endif 
        endif
    endif
enddo

if (A/=0 .and. nAgeClasses>1) then  ! A may be not-yet-assigned Parent of Sib f
    do f=1, nS(SB,kB)
        if (MaybeFS(f)/=1 .or. Par(f)/=0 .or. Parent(SibID(f, SB, kB), 3-kB)/=0) cycle  
        if (AgeDiff(SibID(f, SB, kB), A) <= 0) cycle
        call CalcU(-SB, kB, A, kB, LLtmp(1))
        Parent(SibID(f, SB, kB), 3-kB) = A
        call CalcU(-SB, kB, A, kB, LLtmp(2))
        Parent(SibID(f, SB, kB), 3-kB) = 0
        call CalcCLL(SB, kB)
        if (LLtmp(1) - LLtmp(2) < thLRrel) then
            MaybeFS(f) = 0
        endif
    enddo
endif
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
    LL = 777
    return
endif
dLL = 999
do l=1,nSnp
    do f=1, nS(SB,kB)
        if (MaybeFS(f) /= 1) cycle
        do x=1,3
            PrX(x,:) = XPr(2,x,l, SB, kB)
            do g=1,nS(SB,kB)
                Bj = SibID(g, SB, kB)
                if (NFS(Bj) == 0) cycle
                if (g==f) then
                    call ParProb(l, Par(g), 3-kB, -1,0, PrY(:,1))
                else
                    call ParProb(l, Parent(Bj, 3-kB), 3-kB, Bj,-1, PrY(:,1))
                endif
                PrY(:,2) = PrY(:,1)  ! 1: FS, 2: not FS (incl. HS via par(f))
                do y=1,3 
                    do i=1,nFS(Bj)
                        if (Genos(l,FSID(i, Bj))/=-9) then
                            PrY(y,:) = PrY(y,:) * OKA2P(Genos(l,FSID(i,Bj)), x, y, l)
                        endif
                    enddo
                    if (g==f) then
                        if (Par(g) < 0) then 
                            do h = 1, nS(-Par(g), 3-kB)
                                Ei = SibID(h, -Par(g), 3-kB)                               
                                if (Parent(Ei, kB) == -SB) cycle  
                                if (NFS(Ei) == 0) cycle  
                                call ParProb(l, Parent(Ei,kB), kB, Ei,-1, PrZ)
                                do z=1,3                 
                                    do i=1, nFS(Ei)
                                        if (FSID(i,Ei) == A) cycle
                                        if (Genos(l, FSID(i,Ei))/=-9) then
                                            PrZ(z) = PrZ(z) * OKA2P(Genos(l,FSID(i,Ei)),y,z,l)
                                        endif
                                    enddo
                                enddo
                                PrY(y,:) = PrY(y,:) * SUM(PrZ)
                            enddo
                        endif
                        if (A/=0) then
                            if (Genos(l, A)/=-9) then
                                PrY(y,1) = PrY(y,1) * OKA2P(Genos(l,A), x, y, l)
                                if (PA/=0) then
                                    PrY(y,2) = PrY(y,2) * OKAP(Genos(l,A), y, l)
                                else
                                    PrY(y,2) = PrY(y,2) * OHWE(Genos(l,A), l)
                                endif
                            endif
                        else if (SA/=0) then
                            PrY(y,1) = PrY(y,1) * SUM(XPr(1,:,l, SA,kA) * AKA2P(:, x, y))
                            if (PA/=0) then
                                PrY(y,2) = PrY(y,2) * SUM(XPr(1,:,l, SA,kA) * AKAP(:, y, l))
                            else
                                PrY(y,2) = PrY(y,2) * SUM(XPr(1,:,l, SA,kA) * AHWE(:,l))
                            endif
                        endif 
                    endif
                enddo  ! y
                do i=1,2
                    PrX(x,i) = PrX(x,i) * SUM(PrY(:,i))
                enddo
            enddo  ! g
        enddo  ! x
        PrL(l,f,:) = LOG10(SUM(PrX, DIM=1))
    enddo  ! f   
enddo

dLL = 777
do f = 1, nS(SB, kB)
    if (MaybeFS(f)/=1) cycle
    dLL(f) = SUM(PrL(:, f,1)) - SUM(PrL(:,f,2))
enddo

if (A/=0) then
    call CalcU(A,kA, -SB, kB, LLUX)
    LL = MAXVAL(dLL, MASK=dLL/=777) + LLUX   
else if (SA/=0) then
    call CalcU(-SA, kA, -SB, kB, LLUX)
    do f = 1, nS(SB, kB)
        if (MaybeFS(f)/=1) cycle
        if (Par(f)==0) then
            dLL(f) = SUM(PrL(:, f, 1)) - SUM(PrL(:, f, 2)) + LLUX
        else  ! consider changes in SA (e.g. inbreeding loops) 
            GpID(3-kB, SA, kA) = Par(f)  ! temporary assigned shared parent w. BB(f)
            call PairUA(-SA, -SB, kA, kB, dLL(f))  ! additionally parent via SB
        endif
    enddo
    GpID(3-kB, SA, kA) = PA
    LL = MaxLL(dLL) 
endif
    
end subroutine AddFS

! ##############################################################################################

subroutine AddParent(A, SB, k, LL)  ! is A parent of sibship SB?
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l,x,y,m, G(2)
double precision :: PrL(nSnp), PrX(3), PrXY(3,3), PrG(3, 2)

LL = 999
G = 0
if (ANY(AgeDiff(SibID(1:nS(SB,k),SB,k), A)<=0)) then  ! A too young to be P or GP of i
    LL = 777
    return
endif

do m=1,2
    if (Parent(A,m)/= 0) then   ! todo: allow for sibship/real parent
        if (GpID(m,SB,k)/= 0 .and. GpID(m,SB,k) /= Parent(A,m)) then
            LL = 777
            return
        else
            G(m) = Parent(A,m)
        endif
    else if(GpID(m,SB,k)/=0) then
        G(m) = GpID(m,SB,k)
    endif
enddo

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9) then
        PrL(l) = LOG10(SUM(XPr(3,:,l, SB,k)))
    else
        call ParProb(l, A, k, 0,0, PrX)
        do m=1,2
            call ParProb(l, G(m), m, 0,0, PrG(:,m))    
        enddo
        do x=1,3
            do y=1,3
                PrXY(x,y) = XPr(1,x,l, SB,k) * PrX(x) * &
                  SUM(AKA2P(x, y, :) *  PrG(y,3-k) * PrG(:, k))
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXY))
    endif
enddo

LL = SUM(PrL)

end subroutine AddParent

! ##############################################################################################

subroutine AddGP(A, SB, k, LL)  ! add A as a grandparent to sibship SB
! TODO: check if sharing a parent 3-k
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, m, i, cat, curGP
double precision :: PrL(nSnp), PrY(3), PrX(3), LLtmp(3)

LL = 999
if (Sex(A)/=3) then
    m = Sex(A)
!    if (GpID(m,SB,k) /= 0) then  ! allow replacement  (else change CalcGPZ)
!        LL = 777
!    endif
else if (GpID(1,SB,k)==0) then
    m = 1
else if (GpID(2,SB,k)==0) then
    m = 2
else
    LL = 777
endif
if (LL==777) return

 cat = 0
if (GpID(3-k,SB,k) < 0) then
    if (Parent(A, 3-k)==GpID(3-k,SB,k)) then
        cat = 1
    else
        do i=1,nS(SB,k)
            if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then 
                cat = 1
                exit
            endif
        enddo
    endif
endif

do i=1, nS(SB, k)
    if (AgeDiff(A, SibID(i,SB,k))==999) cycle
    if (AgeDiff(SibID(i,SB,k), A)<=0) then  ! A too young to be P or GP of i
        LL = 777
    else if (k==1 .and. m==1 .and. AgePriorM(AgeDiff(SibID(i,SB,k), A)+1, 3)==0.0) then
        LL = 777
    else if (k==2 .and. m==2 .and. AgePriorM(AgeDiff(SibID(i,SB,k), A)+1, 4)==0.0) then
        LL = 777 
    else if (k/=m .and. AgePriorM(AgeDiff(SibID(i,SB,k), A)+1, 5)==0.0) then
        LL = 777         
    endif
enddo
if (LL==777) return

if (cat==0) then
    PrL = 0
    do l=1,nSnp
        if (Genos(l,A)==-9) then
            PrL(l) = LOG10(SUM(XPr(3,:,l, SB,k)))
        else
            call ParProb(l, GpID(3-m, SB, k), 3-m, 0,0, PrY)
            do x=1,3
                PrX(x) = XPr(1,x,l, SB,k) * SUM(AKOAP(x, Genos(l,A), :, l) * PrY)
            enddo
        endif
        PrL(l) = LOG10(SUM(PrX))
    enddo
    LL = SUM(PrL) + Lind(A)
else  ! inbreeding loop present / will be created
    if (GpID(3-m, SB, k) < 0) then
        call CalcU(-SB, k, GpID(3-m, SB, k), 3-m, LLtmp(1))
    endif
    curGP = GPID(m, SB, k)
    GpID(m, SB, k) = A
    call CalcU(-SB, k, A, 3-k, LL)
    if (GpID(3-m, SB, k) < 0) then
        call CalcU(-SB, k, GpID(3-m, SB, k), 3-m, LLtmp(2))
        LL = LL + (LLtmp(2) - LLtmp(1))
    endif
    GPID(m,SB,k) = CurGP
    call CalcCLL(SB, k)
endif
      
end subroutine AddGP

! ##############################################################################################

subroutine AddGGP(A, SB, k, LL)
use Global
implicit none
! A a great-grandparent of sibship SB? (only calculating over non-grandparent-assigned)

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, m, y, AncG(2,16), i
double precision :: PrL(nSnp), PrXY(3,3), PrZ(3)

LL = 999
if (GpID(1, SB,k)/=0) then
    if (GpID(2, SB,k)/=0) then  ! should be assigned as parent-of-gp
        LL = 777   !(or 888)
    else
        m = 2
    endif
else
    m = 1  ! doesn't really matter (?); GpID(m, SB, k) == 0.
endif
if (LL==777) return

if (GpID(3-m, SB, k)/=0) then
    call GetAncest(GpID(3-m, SB, k), 3-m, AncG)   
    if (ANY(AncG == A)) then
        LL = 444  ! possible; not yet implemented
    else if ((Parent(A,1)/=0 .and. ANY(AncG(1, 2:16) == Parent(A,1))) .or. &
      (Parent(A,2)/=0 .and. ANY(AncG(2, 2:16) == Parent(A,2)))) then
        LL = 444  ! possible; not yet implemented
    endif
endif
if (GpID(3-k,SB,k) < 0) then
    do i=1,nS(SB,k)
        if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then 
            LL = 444
            exit
        endif
    enddo
endif
if (LL==444) return ! TODO: implement.
PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9) then
        PrL(l) = LOG10(SUM(XPr(3,:,l, SB,k)))
    else
        call ParProb(l, GpID(3-m, SB, k), 3-m, 0,0, PrZ)
        do x=1,3  ! sibship parent
            do y=1,3
                PrXY(x,y) = XPr(1,x,l, SB,k) * SUM(AKA2P(x, y, :) * PrZ &
                * AKOP(y, Genos(l,A), l))
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXY))
    endif            
enddo
LL = SUM(PrL) + Lind(A)
     
end subroutine AddGGP

! ##############################################################################################

subroutine ParentHFS(A, SA, kA, SB, kB, hf, LL)  ! parents of SA and SB HS/FS?
use Global
implicit none

integer, intent(IN) :: A, SA, kA, SB, kB, hf
double precision, intent(OUT) :: LL
integer :: m, G(2), l, x, y, u,v, AncA(2,16), AncB(2,16), i, j,z, r, Ei, GA, GB,e, &
    DoneA(MaxSibSize), Ai, Bj, nA, AA(maxSibSize), catA(maxSibSize), catB(nS(SB,kB)+1)
double precision :: PrG(3,2), PrL(nSnp), PrXV(3,3,3,3,3,2), PrPA(3, 2), LLm(2), &
    PrGA(3), PrGB(3), PrE(3,2), PrH(3)

LLm = 999
G = 0  

 if (A/=0) then
    call GetAncest(A, kA, AncA)
    nA = 1
    AA(1) = A
else
    call GetAncest(-SA, kA, AncA)
    nA = nS(SA,kA)
    AA(1:nA) = SibID(1:nA, SA, kA)
endif
 call GetAncest(-SB, kB, AncB)
 
  G = 0 
do m=1,2
    if (m/=hf .and. hf/=3) cycle
    if (AncA(m, kA+2)/=0) then
        if (AncA(m, kA+2) == -SB) then
            LLm(m) = 777
        else if (AncB(m, kB+2)/=0 .and. AncA(m, kA+2)/=AncB(m, kB+2)) then
            LLm(m) = 777
        else if (AncB(m, kB+2)==0) then
            G(m) = AncA(m, kA+2)
        else if (AncB(m, kB+2)==AncA(m, kA+2)) then
            G(m) = AncA(m, kA+2)
!            LLm(m) = 888  ! already are sibs
        else
            LLm(m) = 777
        endif
    else
        if (AncB(m,kB+2)/=0 .and. AncB(m,kB+2) == AncA(m, 2)) then
            LLm(m) = 777
        else 
            G(m) = AncB(m, kB+2)
        endif
    endif
    if (hf==3) then  ! FS
        if (ANY(AncA(kB, 3:16) == -SB)) then
            LLm = 777
        else if (A>0) then
            if (ANY(AncB == A)) then
                LLm = 777
            endif
        else if (SA/=0) then
            if (ANY(AncB(kA,3:16) == -SA)) then
                LLm = 777
            endif
        endif
    endif
    do x=2,16
        if (AncB(m,x) > 0) then
            if (A > 0) then
                if (AgeDiff(A, AncB(m,x)) < 0) then
                    LLm(m) = 777  ! A older than putative ancestor
                endif 
            else 
                if (ANY(AgeDiff(SibID(1:nS(SA,kA),SA,kA), AncB(m,x)) < 0)) then
                    LLm(m) = 777  ! Ai older than putative ancestor
                endif
            endif
        endif
        if (AncA(m,x) > 0) then
            if (ANY(AgeDiff(SibID(1:nS(SB,kB),SB,kB), AncA(m,x)) < 0)) then
                LLm(m) = 777 
            endif
        endif
    enddo
enddo

if (hf==3) then
    if (LLm(1)==777 .or. LLm(2)==777) then
        LL = 777
    endif
else
    if (LLm(hf)==777) then 
        LL = 777
    else 
        GA = AncA(3-hf, kA+2)
        GB = AncB(3-hf, kB+2)
    endif
endif

if (ANY(AncA(kB, 5:16)==-SB)) then
    LL = 444  ! highly unlikely (but not strictly impossible: TODO)
else if (AncA(kA,2)/=0 .and. ANY(AncB(kA, 5:16) == AncA(kA,2))) then  ! PA
    LL = 777
endif
if (LL /=999) return  
   
 catA = 0  
 catB = 0
 
do i=1, nA
    if (kA/=kB) then
        if (Parent(AA(i), kB) == AncB(kB, 2) .and. AncB(kB, 2)<0) then
            catA(i) = 1
        endif
    else if (kA == kB .and. Parent(AA(i), 3-kA)<0) then  
        do j=1, nS(SB, kB)
            if (Parent(AA(i), 3-kA) == Parent(SibID(j,SB,kB), 3-kB)) then
                catA(i) = 2
                catB(j) = 2
            endif
        enddo
    endif
    if (Parent(AA(i), 3-kA) < 0) then
        if (G(3-kA) == Parent(AA(i), 3-kA)) then  ! incl. hf==3
            if (kA==kB) then
                catA(i) = 3  ! (u) 3-kA = 3-kB == hf 
            else if (kA/=kB) then
                catA(i) = 4  ! (z)
            endif 
        else if (kA==hf .and. GA == Parent(AA(i), 3-kA)) then
            catA(i) = 4  ! (z)
        else if (kA==hf .and. GB == Parent(AA(i), 3-kA)) then
            catA(i) = 5  ! (v)
        endif
    endif
enddo    

do i=1, nS(SB, kB)
    if (kA/=kB) then
        if (Parent(SibID(i,SB,kB), kA) == AncA(kA, 2) .and. AncA(kA,2)<0) then
            catB(i) = 1
        endif
    endif
    if (Parent(SibID(i,SB,kB), 3-kB) < 0) then
        if (G(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
            catB(i) = 3  ! (u)  (for hf<3 .and. hf==3)
        else if (kB==hf .and. GA == Parent(SibID(i,SB,kB), 3-kB)) then
            catB(i) = 4  ! (z) (GA of type 3-kB if hf==kB) (hf==3: G_kB cannot be opp. parent)
        else if (kB==hf .and. GB == Parent(SibID(i,SB,kB), 3-kB)) then
            catB(i) = 5  ! (v)
        endif
    endif
enddo 

PrL = 0
do l=1,nSnp
    do m=1,2
        if (m/=hf .and. hf/=3) cycle
        if (ANY(CatA==3) .or. ANY(CatB==3)) then
            call ParProb(l, G(m), m, -1,0, PrG(:,m)) 
        else
            call ParProb(l, G(m), m, 0,0, PrG(:,m)) 
        endif
    enddo
    if (hf < 3) then
        if (ANY(CatA==4) .or. ANY(CatB==4)) then
            call ParProb(l, GA, 3-hf, -1,0, PrGA)
        else
            call ParProb(l, GA, 3-hf, 0,0, PrGA)
        endif
        if (ANY(CatA==5) .or. ANY(CatB==5)) then
            call ParProb(l, GB, 3-hf, -1,0, PrGB)
        else
            call ParProb(l, GB, 3-hf, 0,0, PrGB)
        endif
    endif
    if (A>0) then
        if (Genos(l,A)==-9) then
            PrL(l) = LOG10(SUM(XPr(3,:,l, SB,kB)))
            cycle
        else  ! TODO: PrPA for Parent(A,kA) if /=0?
            do m=1,2
                call ParProb(l, Parent(A,m), m, A,0, PrPA(:,m))
            enddo
        endif
    endif
    
    PrXV = 0
    do x=1,3  ! SA/PA
        do y=1,3  ! SB
            do u=1,3  ! G_hf / G_3-kB (hf==3)
                do z=1,3  ! G_A (hf/=3) / G_kB (hf==3)
                    do v=1,3 ! G_B (hf/=3)
                        if (hf==3) then
                            PrXV(x,y,u,z,z,:) = AKA2P(x,u,z) * AKA2P(y,u,z) * PrG(u,kB) * PrG(z, 3-kB)
                        else
                            if (GA < 0 .and. GA == -SB) then
                                PrXV(x,y,u,y,v,:) = AKA2P(x,u,y) * AKA2P(y,u,v) * PrG(u,hf) * PrGB(v)
                            else if (GB < 0 .and. GB == -SA) then
                                PrXV(x,y,u,z,x,:) = AKA2P(x,u,z) * AKA2P(y,u,x) * PrG(u,hf) * PrGA(z)
                            else
                                PrXV(x,y,u,z,v,:) = AKA2P(x,u,z) * AKA2P(y,u,v) *PrG(u,hf) * PrGA(z) * PrGB(v)
                            endif
                        endif
                        if (A /=0) then
                            if (Parent(A, kA)/=0) then
                                PrXV(x,y,u,z,v,:) = PrXV(x,y,u,z,v,:) * PrPA(x, kA)
                            endif
                        endif
                    enddo
                enddo
            enddo
            
            DoneA = 0            
            if (ALL(catA==0) .and. ALL(catB==0)) then
                if (SA/=0) then
                    PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * XPr(1,x,l, SA,kA) * XPr(1,y,l, SB,kB)
                else if (A>0) then
                    PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(OKA2P(Genos(l,A), x, :, l) * PrPA(:,3-kA)) * &
                      XPr(1,y,l, SB,kB)
                endif               
        
            else
                do r=1, nS(SB,kB)
                    Bj = SibID(r, SB, kB) 
                    if (NFS(Bj) == 0) cycle 
                    if (catB(r)==0 .or. catB(r)==2) then
                        call ParProb(l, Parent(Bj,3-kB), 3-kB, -1,0, PrE(:,1))
                        PrE(:,2) = PrE(:,1)  ! 1=all; 2=non-sibs only
                    else
                        PrE = 1
                    endif                    
         
                    do j=1, nFS(Bj)  
                        if (Genos(l,FSID(j, Bj))/=-9) then
                            PrE(:,1) =  PrE(:,1) * OKA2P(Genos(l,FSID(j,Bj)), y, :, l)
                        endif
                    enddo

                    if (catB(r)==2) then  ! kA==kB, share parent 3-kB
                        do v = 1, nA
                            Ai = AA(v)
                            if (SA/=0 .and. nFS(Ai) == 0) cycle
                            if (Parent(Ai, 3-kA)/=Parent(Bj,3-kB)) cycle
                            do i=1, nFS(Ai)
                                if (A/=0 .and. FSID(i, Ai)/=A) cycle
                                if (Genos(l,FSID(i, Ai))==-9) cycle
                                PrE(:,1) =  PrE(:,1) * OKA2P(Genos(l,FSID(i,Ai)), x, :, l)
                            enddo
                            doneA(v) = 1
                        enddo
                    endif
                        
                    if (Parent(Bj,3-kB) <0 .and. CatB(r)/=1) then
                        do v=1, nS(-Parent(Bj,3-kB), 3-kB)
                            Ei = SibID(v, -Parent(Bj,3-kB), 3-kB)
                            if (nFS(Ei) == 0) cycle
                            if (Parent(Ei, kB) == -SB) cycle
                            if (Parent(Ei, kA) == Parent(AA(1),kA) .and. Parent(AA(1),kA)/=0) cycle
                            if (A>0 .and. Ei==Parent(AA(1),kA)) cycle
                            do e=1,3
!                                if (Parent(Ei, kB)==G(kB) .and. G(kB)/=0) then  
!                                    PrH = 1
!                                else
                                    call ParProb(l, Parent(Ei, kB), kB, Ei,-1, PrH) ! TODO: Parent(Ei, kB)==GPs
 !                               endif
                                do i=1, nFS(Ei)
                                    if (Genos(l,FSID(i, Ei))==-9) cycle
                                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e, l)
                                enddo
                                PrE(e,:) = PrE(e,:) * SUM(PrH)
                            enddo  
                        enddo
                    endif
                    
                    do i=1,2
                        if (catB(r)==0 .or. catB(r)==2) then 
                            PrXV(x,y,:,:,:,i) = PrXV(x,y,:,:,:,i) * SUM(PrE(:,i))
                        else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
                            PrXV(x,y,:,:,:,i) = PrXV(x,y,:,:,:,i) * PrE(x,i)
                        else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
                            do u=1,3
                                PrXV(x,y,u,:,:,i) = PrXV(x,y,u,:,:,i) * PrE(u,i)
                            enddo
                        else if (catB(r)==4) then  ! hf/=kB, Parent(Bj,3-kB) = G
                            do z=1,3
                                PrXV(x,y,:,z,:,i) = PrXV(x,y,:,z,:,i) * PrE(z,i)
                            enddo
                        else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
                            do v=1,3
                                PrXV(x,y,:,:,v,i) = PrXV(x,y,:,:,v,i) * PrE(v,i)
                            enddo
                        endif
                    enddo
                enddo
                
                do r = 1, nA
                    if (doneA(r)==1) cycle
                    if (SA/=0 .and. NFS(AA(r)) == 0) cycle
                    if (kA/=kB .and. Parent(AA(r),3-kA)==-SB) cycle  ! done
                    if (catA(r)==0) then
                        call ParProb(l, Parent(AA(r),3-kA), 3-kA, -1,0, PrE(:,1))
                        PrE(:,2) = PrE(:,1)
                    else
                        PrE = 1
                    endif
                    
                    do e=1,3
                        do i=1, nFS(AA(r))
                            if (Genos(l,FSID(i, AA(r)))==-9) cycle
                            if (SA/=0 .or. FSID(i, AA(r))==A) then 
                                PrE(e,1) =  PrE(e,1) * OKA2P(Genos(l,FSID(i,AA(r))), x, e, l)
                                doneA(r) = 2
                            else
                                PrE(e,:) =  PrE(e,:) * OKA2P(Genos(l,FSID(i,AA(r))), x, e, l)
                            endif
                        enddo
                                          
                        if (Parent(AA(r), 3-kA) < 0 .and. (SA/=0 .or. ANY(FSID(1:nFS(AA(r)), AA(r))==A))) then
                            do i=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                                Ei = SibID(i, -Parent(AA(r), 3-kA), 3-kA)
                                if (nFS(Ei) == 0) cycle
                                if (Parent(Ei, kB) == -SB) cycle
                                if (SA/=0 .and. Parent(Ei, kA) == -SA) cycle
                                call ParProb(l, Parent(Ei, kA), kA, Ei,-1, PrH)  ! Assume is not one of GPs
                                do j=1, nFS(Ei)
                                    if (A/=0 .and. FSID(i, Ei)==A) cycle
                                    if (Genos(l,FSID(j, Ei))/=-9) then
                                        PrH = PrH * OKA2P(Genos(l,FSID(j,Ei)), :, e, l)
                                    endif
                                enddo
                                PrE(e,:) = PrE(e,:) * SUM(PrH)
                            enddo
                        endif
                    enddo
                    
                    do i=1,2
                        if (catA(r)<3) then
                            PrXV(x,y,:,:,:,i) = PrXV(x,y,:,:,:,i) * SUM(PrE(:,i))
                         else if (catA(r)==3) then 
                            do u=1,3
                                PrXV(x,y,u,:,:,i) = PrXV(x,y,u,:,:,i) * PrE(u,i)
                            enddo
                        else if (catA(r)==4) then 
                            do z=1,3
                                PrXV(x,y,:,z,:,i) = PrXV(x,y,:,z,:,i) * PrE(z,i)
                            enddo
                        else if (catA(r)==5) then 
                            do v=1,3
                                PrXV(x,y,:,:,v,i) = PrXV(x,y,:,:,v,i) * PrE(v,i)
                            enddo
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo

    PrL(l) = LOG10(SUM(PrXV(:,:,:,:,:,1)) / SUM(PrXV(:,:,:,:,:,2)))
enddo

LL = SUM(PrL)

end subroutine ParentHFS

! ##############################################################################################

subroutine DummyGP(SA, SB, kA, kB, LL)  ! SB GP of SA? via observed or unobserved B
! TODO: currently assumes SA and SB are independent  -> too low when sharing opp. sibship
use Global
implicit none

integer, intent(IN) :: SA, SB, kA, kB
double precision, intent(OUT) :: LL
integer :: i, m, l, x, y, z, InBoth(nS(SB,kB)), G(2), v, w, Bi, r, AncB(2,16)
double precision :: LLGX(2), PrL(nSnp), PrZ(3), PrXYZ(3,3,3,3), PrG(3), &
    PrPG(3), PrYV, PrW(3), LLtmp(maxSibSize, 2)

InBoth = 0
do i=1, nS(SB,kB)
    if (kA /= kB) then
        if (Parent(SibID(i,SB,kB), 3-kB) == -SA) then
            LL = 444
        endif
    else if (kA == kB) then
        do r= 1, nS(SA, kA)
            if (Parent(SibID(i,SB,kB), 3-kB)==Parent(SibID(r,SA,kA), 3-kA) .and. &
              Parent(SibID(i,SB,kB), 3-kB)/=0) then
                LL = 444
            endif
        enddo
    endif
enddo
if (LL == 444)  return  ! TODO

 call GetAncest(-SB, kB, AncB)
if (ANY(AncB(kA,3:16) == -SA)) then
    LL = 777 
    return
endif
G = GpID(:, SA, kA)

LLGX = 999
LLtmp = 999
do m=1,2
    if (m==kB .and. GpID(kB, SA, kA) == -SB) then
        LLGX(m) = 777
        cycle
    else if (G(m) > 0) then
        call AddSib(G(m), SB, kB, LLtmp(1,m))
        call AddFS(G(m), SB, kB,0,m, LLtmp(2,m))
        if (MaxLL(LLtmp(:,m)) < 0) then
            LLGX(m) = MaxLL(LLtmp(:,m)) + CLL(SA, kA)
            if (Parent(G(m), kB) /= -SB) then
                LLGX(m) = LLGX(m) - Lind(G(m))
            endif
        else
            LLGX(m) = MaxLL(LLtmp(:,m))
        endif
        cycle
    else if (G(m) == 0) then
        do r=1, nS(sB,kB)
            Bi = SibID(r, sB, kB) 
            if (Sex(Bi)/=m .and. Sex(Bi)/=3) cycle
            call AddGP(Bi, SA, kA, LLtmp(r,m))
            if (LLtmp(r,m) < 0) then
                LLtmp(r,m) = LLtmp(r,m) - Lind(Bi) + CLL(SB,kB)
            endif
        enddo
        LLGX(m) = MaxLL(LLtmp(:,m))
    endif
    PrL = 0
    do l=1,nSnp
         call ParProb(l, GpID(3-m, SA, kA), 3-m, 0,0, PrZ)
        if (G(m)/=0) then
            call ParProb(l, G(m), m, 0,0, PrG)
            if (G(m) > 0) then
                call ParProb(l, Parent(G(m), 3-kB), 3-kB, 0,0, PrPG)
            else if (G(m) < 0) then
                call ParProb(l, GpID(3-kB, -G(m), m), 3-kB, 0,0, PrPG)
            endif
        endif
        PrXYZ = 0
        do x=1,3  ! SA
            PrYV = 1
            do y=1,3  ! parent of SA, offspr of SB
                if (G(m)==0) then 
                    do z=1,3  ! other parent of SA
                        do v=1,3  ! SB
                            PrXYZ(x,y,z,v) = XPr(1,x,l, SA,kA) * AKA2P(x, y, z) * &
                            PrZ(z) * AKAP(y, v, l) * XPr(3,v,l, SB,kB)
                        enddo
                    enddo
                else  
                    do v=1,3  ! SB
                        do r=1, nS(sB,kB)
                            Bi = SibID(r, sB, kB)  
                            if (NFS(Bi) == 0) cycle                           
                            if (ANY(FSID(1:nFS(Bi), Bi) == G(m))) then
                                PrW = PrPG
                            else
                                call ParProb(l, Parent(Bi, 3-kB), 3-kB, Bi,-1, PrW)
                            endif
                            do w=1,3
                                do i=1, nFS(Bi)
                                    if (FSID(i,Bi) == G(m)) cycle
                                    if (Genos(l,FSID(i, Bi))/=-9) then
                                        if (FSID(i,Bi) == G(m)) then
                                            PrW(w)  = PrW(w) * PrG(y) * AKA2P(y, v, w)
                                        else
                                            PrW(w) = PrW(w) * OKA2P(Genos(l,FSID(i,Bi)), v, w, l)
                                        endif
                                    endif
                                enddo
                            enddo
                            PrYV = PrYV * SUM(PrW)
                        enddo
                        do z=1,3
                        PrXYZ(x,y,z,v) = XPr(1,x,l, SA,kA) * AKA2P(x, y, z) * &
                            PrYV * PrZ(z) *  XPr(2,v,l, SB, kB)
                        enddo
                    enddo  
                endif
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXYZ))         
    enddo
    LLGX(m) = MaxLL((/ SUM(PrL),  LLGX(m)/))
enddo

LL = MaxLL(LLGX)
end subroutine DummyGP
! ##############################################################################################

subroutine CalcAgeLR(s, k, CandSib, CandGP, CandCl, CandUA, ALR)
use Global
implicit none

integer, intent(IN) :: s, k, CandSib, CandGP, CandCl, CandUA
double precision, intent(OUT) :: ALR
integer :: i, j

ALR = 0.0
do i = 1, nS(s,k)
    if (BY(SibID(i,s,k))<0) cycle
    if (CandSib /= 0) then
        if (AgeDiff( SibID(i,s,k), CandSib)/=999) then
            ALR = ALR + LOG10(AgePriorM(ABS(AgeDiff( SibID(i,s,k), CandSib))+1, k)) 
        endif
    endif
    if (CandGP /= 0) then
        if (AgeDiff(SibID(i,s,k), CandGP)>=0 .and. AgeDiff(SibID(i,s,k), CandGP)/=999) then
            if (Sex(CandGP)==k) then
                ALR = ALR + LOG10(AgePriorM(AgeDiff(SibID(i,s,k), CandGP)+1, 2+k)) 
            else
                ALR = ALR + LOG10(AgePriorM(AgeDiff(SibID(i,s,k), CandGP)+1, 5)) 
            endif
        endif
    endif
    if (CandCl /= 0) then
        do j=1, nS(CandCl, k)
            if (AgeDiff( SibID(i,s,k), SibID(j,CandCl,k))/=999) then
                ALR = ALR + LOG10(AgePriorM(ABS(AgeDiff( SibID(i,s,k), SibID(j,CandCl,k) ))+1, k)) 
            endif
        enddo
    endif
    if (CandUA /= 0) then
        if (AgeDiff(CandUA, SibID(i,s,k))/=999) then
            ALR = ALR + LOG10(AgePriorM(ABS(AgeDiff(CandUA, SibID(i,s,k)))+1, 6)) 
        endif
    endif
enddo

if (ALR < -HUGE(0.0D0)) then
    ALR = 777
endif

end subroutine CalcAgeLR

! ##############################################################################################

subroutine CalcAgeLRCAU(A, kA, B, kB, ALR)  
! SB old enough to be grandparent of A's? (i.e., est BY(SA) > est. BY(SB) ?)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
double precision, intent(OUT) :: ALR
integer :: i, y, x, BYAB(maxSibSize, 2), mx, nA, nB
double precision :: EstBY(nAgeClasses, 2), ALRtmp(nAgeClasses, nAgeClasses, 2)
! age prior SB parent of SA, using priors on their est. BY, based on their offspring's BY

if (nAgeClasses == 1) then
    ALR = 777
    return
endif

if (A<0) then
    nA = nS(-A,kA)
else if (A>0) then
    nA = 1
endif

if (B<0) then
    nB = nS(-B,kB)
else if (B>0) then
    nB = 1
endif

BYAB = -1
if (A<0) then
    do i = 1, nA
        BYAB(i, 1) = BY(SibID(i,-A,kA))
    enddo
else if (A>0) then
    BYAB(1, 1) = BY(A)    
endif
if (B < 0) then
    do i = 1, nB
        BYAB(i, 2) = BY(SibID(i,-B,kB))
    enddo
else if (B > 0) then
    BYAB(1,2) = BY(B)
endif

mx = COUNT(AgePriorM(:, 6+kB) > 0)  ! max age of parent
if (MINVAL(BYAB, MASK=BYAB>=0) < mx) then
    WHERE(BYAB >= 0) BYAB = BYAB + mx
endif

EstBY = 1
do y=1,nAgeClasses
    do i = 1, nA
        if (BYAB(i,1) < 0) cycle  ! unknown BY
        if (BYAB(i,1) - y < 0 .or. (BYAB(i,1) - y +1) > nAgeClasses) then
            EstBY(y,1) = 0  ! Sib i born prior to year y
            exit
        else
            EstBY(y,1) = EstBY(y,1) * AgePriorM(BYAB(i,1) - y +1, 6+kA)
        endif
    enddo
    
    do i = 1, nB
        if (BYAB(i,2) < 0) cycle  ! unknown BY
        if (BYAB(i,2) - y < 0 .or. (BYAB(i,2) - y +1) > nAgeClasses) then
            EstBY(y,2) = 0  ! Sib i born prior to year y
            exit
        else
            EstBY(y,2) = EstBY(y,2) * AgePriorM(BYAB(i,2) - y +1, 6+kB)
        endif
    enddo
enddo

ALRtmp = 0
do y=1,nAgeClasses-1  ! SB
    do x=1, nAgeClasses  ! SA 
        if (x > y) then
            ALRtmp(x,y, 1) = AgePriorM(x - y +1, 6+kB) * EstBY(x,1) * EstBY(y,2) 
        endif
        ALRtmp(x,y, 2) = EstBY(x,1) * EstBY(y,2)   ! SA & SB unrelated
    enddo
enddo
ALR = LOG10(SUM(ALRtmp(:,:,1))) - LOG10(SUM(ALRtmp(:,:,2)))
 
if (SUM(ALRtmp(:,:,1)) == 0 .or. SUM(ALRtmp(:,:,2)) == 0) then
    ALR = 777
endif

end subroutine CalcAgeLRCAU

! ########################################################################################

! ########################################################################################

subroutine BestRel(LLIN, focal, X, dLL)
use Global
implicit none
! return which relationship is most likely
! assuming order PO,FS,HS,GG,FAU,HAU,U in LL vector 
! should exceed others by thLRrel, and exceed U by thLR

double precision, intent(IN) :: LLIN(7)
integer, intent(IN) :: focal
integer, intent(OUT) :: X
double precision, intent(OUT) :: dLL   ! diff best vs next best
double precision :: LL(7), LLtmp(7), LLmax
integer :: i,j, maybe(6)

dLL = 0
maybe = 0
X=0
LL = LLIN

if (focal/=2 .and. LL(2) < 0 .and. LL(2)>=LL(3)) then  ! FS are also PHS/MHS, conditional on 1 parent
    LL(3) = 999
else if (focal==3 .and. LL(3)>LL(2) .and. LL(3)<0) then
    LL(2) = 999   ! want sib vs non-sib
endif

do i=1,6
    if (LL(i)>0) cycle  ! maybe=0
    if ((LL(i) - LL(7)) > thLR) then
        maybe(i) = 1
        do j=1,6
            if (i==j) cycle
            if (LL(j)>0) cycle
            if ((LL(i) - LL(j)) < thLRrel) then
                maybe(i) = 0   ! i has no longer highest LL
            endif
        enddo
    endif
enddo
        
if (SUM(maybe)==0) then
    if ((LL(7) - MAXVAL(LL(1:6), MASK=LL(1:6)<0)) > thLR) then  
        X = 7   ! unrelated
    else
        X = 8   ! unclear
    endif
else if (SUM(maybe)==1) then
    X = MAXLOC(LL(1:6), MASK=maybe==1, DIM=1)
!else  
!    print *, "!!!!!!!"
!    print *, "check BestRel, >1 maybe"
endif

LLtmp = LL
if (MAXLOC(LLtmp, MASK=LLtmp<0, DIM=1) /= 0) then
    LLmax = MAXVAL(LL, MASK=LLtmp<0)
    LLtmp(MAXLOC(LLtmp, MASK=LLtmp<0)) = 999
    dLL = LLmax - MAXVAL(LL, MASK=LLtmp<0)
endif
            
end subroutine BestRel

! ##############################################################################################

subroutine UpdateAllProbs
use Global
implicit none

integer :: i, k, s

do k=1,2
    do s=1,nC(k)
        call CalcCLL(s, k)
    enddo
enddo

do i=1,nInd
    call CalcLind(i)
enddo

end subroutine UpdateAllProbs

! ##############################################################################################

subroutine CalcLind(i)
use Global
implicit none
! assumes CLL is up-to-date

integer, intent(IN) :: i
integer :: l, x, y, k, z
double precision :: PrL(nSnp), Px(3,2), PrXY(3,3),PrXYZ(3,3,3)

PrL = 0
do l=1,nSnp
    do k=1,2
        call ParProb(l, Parent(i,k), k, i,0, Px(:,k))
    enddo
    do x=1,3
        do y=1,3
            if (Genos(l,i)/=-9) then
                PrXY(x,y) = OKA2P(Genos(l,i), x, y,l) * Px(x,1) * Px(y,2)
            else
                do z=1,3
                    PrXYZ(z, x,y) = OKA2P(z, x, y,l) * Px(x,1) * Px(y,2)
                enddo
            endif
        enddo
    enddo
    if (Genos(l,i)/=-9) then
        PrL(l) = LOG10(SUM(PrXY))
    else
        do z=1,3
            LindG(z, l, i) = SUM(PrXYZ(z,:,:))
        enddo
    endif
enddo
                
Lind(i) = SUM(PrL)
LindX(:,i) = PrL

end subroutine CalcLind

! ###########################################################################

subroutine CalcCLL(s, k)
use Global
implicit none
! returns XPr: likelihood;  DumP: probability, scaled  (no age prior.)
! split into 1: sibs only 2: gp effect only, 3: all

integer, intent(IN) :: s, k  ! S: sibship number, k: maternal(1), paternal(2), or unknown(3)
integer :: l, x, i, Ei, r, y, z, g, Ri, v, cat(nS(s,k)+1)
double precision :: PrL(nSnp), PrY(3), PrGG(3,2), PrZ(3), PrXZ(3,3)

 cat = 0
do r=1,nS(s,k)
    Ri = SibID(r, s, k)
    if (Parent(Ri, 3-k)<0 .and. Parent(Ri, 3-k)==GpID(3-k,s,k)) then  ! TODO: /=0 ?
        cat(r) = 1
    endif
enddo

if (ALL(GpID(:,s,k)<0) .and. ALL(cat(1:nS(s,k))==0)) then  ! check if sibship parent is inbred
    if (GPID(1,s,k) == GPID(1, -GPID(2,s,k),2)) then
        cat(nS(s,k)+1) = 2
    else if (GPID(2,s,k) == GPID(2, -GPID(1,s,k),1)) then
        cat(nS(s,k)+1) = 1
    endif
endif

!  nD=0
!do m=1,2
!    do r=1,nC(m)
!        if (GpID(k,r,m) == -s) then
!            nD = nD + 1
!            Ds(nD) = r
!            Dk(nD) = m
!        endif
!    enddo
!enddo

PrL = 0
do l=1,nSnp
    do g=1,2   !grandparents
        if (g/=k .and. ANY(cat(1:nS(s,k))==1)) then
            call ParProb(l, GpID(g,s,k), g, -1,0, PrGG(:,g))
        else if (cat(nS(s,k)+1)==g) then
            PrGG(:,g) = XPr(3,:,l, -GpID(g,s,k),g)
        else
            call ParProb(l, GpID(g,s,k), g, 0,0, PrGG(:,g))
        endif
    enddo
     do x=1,3  ! genotype dummy parent
        do z=1,3
            if (cat(nS(s,k)+1)==0) then
                PrXZ(x,z) = SUM(AKA2P(x,:, z) * PrGG(z,3-k) * PrGG(:,k))  ! GPs
            else if (cat(nS(s,k)+1)==k) then
                PrXZ(x,z) = SUM(AKA2P(x, :, z) * PrGG(:,k))
            else if (cat(nS(s,k)+1)==3-k) then
                PrXZ(x,z) = SUM(AKA2P(x, :, z) * PrGG(z,3-k))
            endif
        enddo
    enddo
    if (cat(nS(s,k)+1)>0) then
        PrXZ = PrXZ/SUM(PrXZ) 
    endif
    
    do x=1,3
        XPr(2,x,l, s,k) = SUM(PrXZ(x,:))  ! GP
        
        do r=1, nS(s,k)
            Ri = SibID(r, s, k)  ! array with IDs
            if (NFS(Ri) == 0) cycle  ! moved to its FS
            if (cat(r)==1) then
                PrY = 1
            else
                call ParProb(l, Parent(Ri, 3-k), 3-k, -1,0, PrY)
            endif
            
            if (Parent(Ri, 3-k) < 0) then
                do y=1,3
                    do v=1, nS(-Parent(Ri, 3-k), 3-k)
                        Ei = SibID(v, -Parent(Ri, 3-k), 3-k)  
                        if (NFS(Ei) == 0) cycle
                        if (Parent(Ei, k) == -s) cycle
                        call ParProb(l, Parent(Ei, k), k, Ei,-1, PrZ)
                        do i=1, nFS(Ei) 
                            if (Genos(l,FSID(i, Ei))==-9) cycle
                            PrZ = PrZ * OKA2P(Genos(l,FSID(i,Ei)), :, y, l)
                        enddo
                        PrY(y) = PrY(y) * SUM(PrZ)
                    enddo
                enddo
                PrY = PrY / SUM(PrY)
            endif
           
            do i=1, nFS(Ri)  ! default: nFS = 1
                if (Genos(l,FSID(i, Ri))==-9) cycle
                do y=1,3
                    PrY(y) = PrY(y) * OKA2P(Genos(l,FSID(i,Ri)), x, y, l)
                enddo
            enddo                
                        
            if (cat(r)==1) then
                PrXZ(x,:) = PrXZ(x,:) * PrY  
            else
                PrXZ(x,:) = PrXZ(x,:) * SUM(PrY)
            endif
        enddo ! r 
        
!        if (nD>0) then  ! Todo: FS dummy - real; GpID(3-k,Ds(r), Dk(r)) == GpID(3-k,s,k)
!            do r=1,nD
!                call ParProb(l, GpID(3-k,Ds(r), Dk(r)), 3-k, 0, PrY(:,2))
!                do y=1,3
!                    PrY(y,2) = PrY(y,2) * SUM(AKA2P(:,x,y) * XPr(1,:,l,DS(r), Dk(r)))
!                enddo
!                do z=1,2
!                    PrXZ(x,:,z) = PrXZ(x,:,z) * SUM(PrY(:,2))
!                enddo
!            enddo
!        endif
    enddo ! x
    do x=1,3  ! account for GP, dumm offspr & connected sibships
        DumP(x,l, s,k) = SUM(PrXZ(x,:))/ SUM(PrXZ)
        XPr(3,x,l, s,k) = SUM(PrXZ(x,:))
        XPr(1,x,l, s,k) = XPr(3,x,l, s,k) / XPr(2,x,l, s,k)  
    enddo 
    PrL(l) = LOG10(SUM(XPr(3,:,l, s,k))) 
enddo
 CLL(s,k) = SUM(PrL)
 
end subroutine CalcCLL

! ##############################################################################################

subroutine ParProb(l, i, k, A, B, prob)
use Global
implicit none

integer, intent(IN) :: l, i, k, A,B
double precision, intent(OUT) :: prob(3)
integer :: x,j, AB(2)
double precision :: PrP(3, 2), PrY(3)

 if (i == 0) then  ! no parent
    prob = AHWE(:, l)
else if (i > 0) then  ! real parent
    if (Genos(l,i) /= -9) then
        prob = AcO(:, Genos(l,i), l)
    else
        prob = LindG(:, l, i)
    endif
else if (i < 0) then  ! dummy parent
    if (A==0) then   ! probability
        prob = DumP(:,l, -i,k)    
    else if (A<0) then  ! grandparent contribution only
        prob = XPr(2,:,l, -i, k)  
    else if (A>0) then   ! exclude indiv A from calc
        if ((Genos(l,A)==-9 .and. (nFS(A)<=1 .or. B>=0)) .or. Parent(A,k)/=i) then
            prob = DumP(:,l, -i,k)
        else
            AB = (/ A, B /)
            do j=1,2
                if (j==2 .and. B<=0) cycle
                if (Parent(AB(j), 3-k)==0) then  
                    PrP(:,j) = AHWE(:,l)
                else if (Parent(AB(j), 3-k)>0) then
                    if (Genos(l,Parent(AB(j), 3-k)) /= -9) then
                        PrP(:,j) = AcO(:, Genos(l, Parent(AB(j),3-k)), l)
                    else
                        PrP(:,j) = LindG(:, l, Parent(AB(j), 3-k))
                    endif
                else if (Parent(AB(j), 3-k)<0) then  
                    PrP(:,j) = DumP(:,l, -Parent(AB(j),3-k), 3-k)  ! TODO: recursive?
                endif
            enddo

            do x = 1, 3
                if (B>=0 .or. nFS(A)<=1) then
                    prob(x) = XPr(3,x,l, -i, k) / SUM(OKA2P(Genos(l,A), x, :, l) * PrP(:,1))    
                else if (B==-1) then  ! exclude all FS of A
                    PrY = PrP(:,1)
                    do j=1, nFS(A)
                        if (Genos(l,FSID(j, A))==-9) cycle
                        PrY = PrY * OKA2P(Genos(l,FSID(j, A)), x, :, l)
                    enddo
                    prob(x) = XPr(3,x,l, -i, k) / SUM(PrY)
                endif
                if (B>0) then
                    if (Genos(l,B)==-9) cycle
                    prob(x) = prob(x) / SUM(OKA2P(Genos(l,B), x, :, l) * PrP(:,2)) 
                endif
            enddo
            prob = prob/SUM(prob)
       endif
    endif
endif

end subroutine ParProb

! ##############################################################################################

subroutine Connected(A, kA, B, kB, Con)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
logical, intent(OUT) :: Con
integer :: i, j, m, nA, nB, AA(maxSibsize), BB(maxSibsize), n

 Con = .FALSE.
if (A==0 .or. B==0) then
    Con = .FALSE.
    return
endif

if (A>0) then
    nA = 1
    AA(1) = A
else
    nA = nS(-A,kA)
    AA(1:nA) = SibID(1:nA, -A, kA)
endif

if (B>0) then
    nB = 1
    BB(1) = B
else
    nB = nS(-B,kB)
    BB(1:nB) = SibID(1:nB, -B, kB)
endif

do j=1, nB
    do i=1, nA
        do m=1,2  
            if (Parent(AA(i), m) < 0) then
                if (Parent(AA(i),m) == Parent(BB(j),m)) then
                    Con = .TRUE.
                    return
                else if(ANY(GpID(:,-Parent(AA(i), m),m) == BB(j))) then
                    Con = .TRUE.
                    return
                else if(ANY(GpID(:,-Parent(AA(i), m),m) < 0)) then
                    do n=1,2
                        if (GpID(n,-Parent(AA(i), m),m) == Parent(BB(j),n) .and. &
                          Parent(BB(j),n)<0) then
                            Con = .TRUE.
                            return
                        endif 
                    enddo
                endif
            endif
            if (Parent(BB(j),m)<0) then
                if (ANY(GpID(:,-Parent(BB(j),m),m) == AA(i))) then
                    Con = .TRUE.
                    return
                else if(ANY(GpID(:,-Parent(BB(j),m),m) < 0)) then
                    do n=1,2
                        if (GpID(n,-Parent(BB(j), m),m) == Parent(AA(i),n) .and. &
                          Parent(AA(i),n)<0) then
                            Con = .TRUE.
                            return
                        endif 
                    enddo
                endif
            endif
        enddo
    enddo
enddo

end subroutine Connected

! ##############################################################################################

subroutine GetAncest(A, k, Anc)
use Global
implicit none

integer, intent(IN) :: A, k
integer, intent(OUT) :: Anc(2, 16)  ! 4 generations back
integer :: m, j, i

 Anc = 0
if (A==0) return

if (A > 0) then  
    if (Sex(A)/=3) then
        Anc(Sex(A),1) = A
    else
        Anc(1, 1) = A
    endif
     Anc(:, 2) = Parent(A, :)
else if (A < 0) then
    Anc(k, 2) = A  
endif

if (ALL(Anc(:, 2)==0)) return
do m = 1, 2 
    if (Anc(m, 2) > 0) then
        Anc(:, 2+m) = Parent(Anc(m, 2), :)
    else if (Anc(m, 2) < 0) then
        Anc(:, 2+m) = GpID(:, -Anc(m, 2), m)
    endif
enddo
if (ALL(Anc(:, 2:3) == 0)) return
do j = 2, 8
    do m = 1, 2
        i = 2 * (j-1) + m
        if (Anc(m, j) > 0) then
            Anc(:, i) = Parent(Anc(m, j), :)
        else if (Anc(m, j) < 0) then
            Anc(:, i) = GpID(:, -Anc(m, j), m)  
        endif
    enddo
enddo

end subroutine GetAncest

! ##############################################################################################

subroutine CalcParentLLR  ! Calc parental LLR (vs next most likely relationship)
use Global
implicit none

integer :: i, k, s, CurPar(2), m, nonG(6), CurGP(2), g
double precision :: LLA(7), LLtmp(2,2,2), LLCP(2,2)

LR_parent = 999 
do i=1, nInd
    if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle   
    CurPar = Parent(i,:)
    Parent(i,:) = 0
    call CalcLind(i)
    do k=1,2  ! remove i from sibgroup
        if (CurPar(k)<0) then
            call RemoveSib(i, -CurPar(k), k)
        endif
    enddo
    
    LLtmp = 999
    do m=1,2  ! m=1: no opp. sex parent;  m=2: with opp. sex parent
        if (m==2 .and. (CurPar(1)==0 .or. CurPar(2)==0)) cycle
        do k=1,2
            if (m==1 .and. CurPar(k) == 0) cycle
            if (m==2) then  ! temp. assign parent 3-k
                call CalcU(i, 3-k, CurPar(3-k), 3-k, LLCP(k,1))
                LLCP(k,1) = LLCP(k,1) - Lind(i)
                Parent(i, 3-k) = CurPar(3-k)
                if (CurPar(3-k)<0) then
                    call DoAdd(i, -CurPar(3-k), 3-k)
                endif
                call CalcU(i, 3-k, CurPar(3-k), 3-k, LLCP(k,2))
                LLCP(k,2) = LLCP(k,2) - Lind(i)
            endif
            
            if (CurPar(k) > 0) then
                call calcPair(i, CurPar(k), k, .FALSE., LLA, 1)
                LLtmp(1,k,m) = LLA(1)
                LLtmp(2,k,m) = MaxLL(LLA(2:7))
            else if (CurPar(k) < 0) then
                if (nS(-CurPar(k),k)==1) then  ! was sibling pair
                    call calcPair(i, SibID(1,-CurPar(k),k), k, .FALSE., LLA, 3)
                else
                    call CheckAdd(i, -CurPar(k), k, LLA, 3)
                endif
                if (m==1) LLA(2) = 333   ! FS does not count here.
                LLtmp(1,k,m) =  MaxLL(LLA(2:3))
                LLtmp(2,k,m) =  MaxLL((/LLA(1), LLA(4:7)/))
            endif
            
            if (m==1) then
                LR_parent(i,k) = LLtmp(1,k,m) - LLtmp(2,k,m)
            else if (m==2) then  
                Parent(i, 3-k) = 0
                call CalcLind(i)
                if (CurPar(3-k)<0) then
                    call RemoveSib(i, -CurPar(3-k), 3-k)
                endif
            endif         
        enddo
    enddo
    if (CurPar(1)/=0 .and. CurPar(2)/=0) then
        LR_parent(i,3) = MINVAL(LLtmp(1,:,2) - MAX(LLtmp(1,:,1), LLtmp(2,:,2)))       
    endif

    Parent(i,:) = CurPar  ! restore
    call CalcLind(i)
    do k=1,2
        if (CurPar(k)<0) then
            call DoAdd(i, -CurPar(k), k)
        endif
    enddo
    call CalcLind(i)
enddo

!parents of dummies (Sibship GPs)
nonG = (/1,2,3,5,6,7/)
do k = 1,2
    do s=1, nC(k)
        CurGP = GpID(:, s, k)
        GpID(:, s, k) = 0
        call CalcCLL(s,k)
        
        LLtmp = 999
        do m=1,2
            if (m==2 .and. (CurGP(1)==0 .or. CurGP(2)==0)) cycle
            do g=1,2
                if (m==1 .and. CurGP(g) == 0) cycle
                if (m==2) then  ! temp. assign GP 3-g
                    GpID(3-g, s, k) = CurGP(3-g)
                    call CalcCLL(s,k)
                endif
                
                if (curGP(g) > 0) then
                    call checkAdd(CurGP(g),s,k, LLA, 4)  ! B=GP + CurGP(m)_7
                else if (curGP(g) < 0) then
                    call checkMerge(s, -CurGP(g), k, g, LLA, 4)  
                    if (m==1) then
                        call PairUA(-s, CurGP(g), k, g, LLA(4))  ! not counting SA FS with a B
                    endif
                endif
                LLtmp(1,g,m) = LLA(4)
                LLtmp(2,g,m) = MaxLL(LLA(nonG))  
                
                if (m==1) then
                    LR_GP(g,s,k) = LLtmp(1,g,m) - LLtmp(2,g,m)
                else if (m==2) then  ! reset to 0
                    GpID(3-g, s, k) = 0
                    call CalcCLL(s,k)
                endif 
            enddo
        enddo
        if (CurGP(1)/=0 .and. CurGP(2)/=0) then
            LR_GP(3,s,k) = MINVAL(LLtmp(1,:,2) - MAX(LLtmp(1,:,1), LLtmp(2,:,2)))       
        endif      

        ! restore
        GpID(:,s,k) = CurGP
        call CalcCLL(s, k)
    enddo
enddo

end subroutine CalcParentLLR

! ##############################################################################################

subroutine RemoveSib(A, s, k)  ! removes individual A from sibship s
use Global
implicit none

integer, intent(IN) :: A, s, k
integer :: u, j, p, curFS, ox, o, h, v

do u=1,nS(s,k)   
    j = SibID(u,s,k)  ! 1st one in FS
    if (nFS(j)==0) cycle
    curFS = 0
    do v=1, nFS(j)
        if (FSID(v, j) == A) then  ! drop FS: move FS group to next lowest ID
            if (nFS(j)>1 .and. curFS==0) then
                p = MINVAL(FSID(1:nFS(j), j), MASK=(FSID(1:nFS(j), j)/=A))  ! note: p==j if A/=j
                curFS = p
                ox = 1
                do o=1, nFS(j)
                    if (FSID(o,j)==A) cycle
                    FSID(ox, p) = FSID(o, j)
                    ox = ox+1
                enddo
                nFS(p) = nFS(j)-1
                nFS(A) = 1
                FSID(1,A) = A
                exit
            endif
        endif
    enddo
enddo

do u=1,nS(s,k)
    if (SibID(u,s,k)==A) then
        do h=u, nS(s, k)-1  ! drop HS
            SibID(h, s, k) = SibID(h+1, s, k)
        enddo
        SibID(nS(s,k), s, k) = 0
    endif
enddo
nS(s,k) = nS(s,k)-1
 call calcCLL(s, k)
 call CalcLind(A)
do u=1,nS(s,k)   ! update LL of connected sibships
    j = SibID(u,s,k)
    if (Parent(j,3-k) < 0) then
        call CalcCLL(-Parent(j,3-k), 3-k)
    endif                    
    call CalcLind(j)
enddo
 call calcCLL(s, k)

Parent(A, k) = 0
call CalcLind(A)
 
end subroutine RemoveSib

! ##############################################################################################

! @@@@   INPUT & PRECALC PROB.   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################

subroutine ReadData
use Global
implicit none

integer :: i,j,k,l, minAgeRepro
integer, allocatable, dimension(:) :: SexTmp, ByTmp
integer, allocatable, dimension(:,:) :: GenosR
 character(len=2000) :: GenoFileName, LifehistFileName, dumC !, AgePriorFileName

! print *, "Reading specs ... "
open (unit=1,file="SequoiaSpecs.txt",status="old")
read (1,*) dumC, GenoFileName      ! 1
read (1,*) dumC, LifehistFileName
read (1,*) dumC, nInd              ! 3
read (1,*) dumC, nIndLH
read (1,*) dumC, nSnp              ! 5 
read (1,*) dumC, Er
read (1,*) dumC, MaxMismatch        ! 7
read (1,*) dumC, thLR           
read (1,*) dumC, thLRrel          ! 9    
read (1,*) dumC, minAgeRepro
read (1,*) dumC, ParentageFileName  ! 11
read (1,*) dumC, AgePriorFileName   ! 12
read (1,*) dumC, nAgeClasses
read (1,*) dumC, maxSibSize         ! 14
read (1,*) dumC, Nrounds
read (1,*) dumC, PedigreeFileName   ! 16
read (1,*) dumC, DumPrefix(1)    
read (1,*) dumC, DumPrefix(2)       ! 18 
 close (1)
 
!=================
 
allocate(GenosR(nInd,nSnp))
allocate(Genos(nSnp, nInd))
Genos = -9
allocate(Id(nInd))
Id = "NA"

open (unit=101,file=trim(GenoFileName),status="old")
!read (101, *)  ! header
do i=1,nInd
!    read (101,*) dumI(1), Id(i), dumI(1:4), GenosR(i,:)
    read (101,*)  Id(i), GenosR(i,:)
    do l=1,nSnp
        if (GenosR(i,l)/=-9) then
            Genos(l,i) = GenosR(i,l)+1
        else
            Genos(l,i) = GenosR(i,l)
        endif
    enddo
enddo
 close (101)
 
!=================
! allele frequencies
allocate(AF(nSNP))
do l=1,nSnp
    AF(l)=float(SUM(GenosR(:,l), MASK=GenosR(:,l)/=-9))/(2*nInd)
enddo 

!=================
allocate(SexTmp(nIndLH))
allocate(NameLH(nIndLH))
allocate(ByTmp(nIndLH))

allocate(BY(nInd))
BY=-999
allocate(Sex(nInd))
Sex = 3
allocate(AgeDiff(nInd,nInd))
AgeDiff=999
 
open(unit=103, file=trim(LifehistFileName), status="old")
read(103, *)
do k=1, nIndLH
    read(103,*) NameLH(k), SexTmp(k), ByTmp(k)
enddo
 close(103)

 ! rearrange lifehistory info to same order as genotype file
do i=1,nInd
    do k=1,nIndLH
        if(Id(i)==NameLH(k)) then
            if (BYtmp(k)>0) then
                BY(i) = BYtmp(k)
            endif
            if (SexTmp(k)==1 .or. SexTmp(k)==2) then
                Sex(i)=SexTmp(k)
            endif
        endif
    enddo
enddo
WHERE (BY /= -999) BY = BY - MINVAL(BY, MASK=BY>0)

!=================
do i=1, nInd
    do j=1, nInd
        if (BY(i)/=-999 .and. BY(j)/=-999) then
            AgeDiff(i,j) = BY(i) - BY(j)   ! if >0, then j older than i
        endif   
    enddo
enddo

!=================
allocate(AgePriorM(nAgeClasses, 8))
open(unit=102, file=trim(AgePriorFileName), status="old")
read(102,*)   ! header
do k=1,nAgeClasses
    read(102,*) AgePriorM(k,:)
enddo 
 close(102)
 
!=================
! allocate arrays
allocate(Lind(nInd))
Lind = 0
allocate(LindX(nSnp, nInd))
LindX = 0
allocate(FSID(MaxSibSize, nInd))
FSID(1, :) = (/ (i, i=1, nInd) /)
allocate(NFS(nInd))
NFS = 1
allocate(DumP(3,nSnp, nInd/2,2))
DumP = 0
allocate(XPr(3,3,nSNP, nInd/2,2))
XPr = 0  
allocate(GpID(2, nInd/2,2))
GpID = 0 
allocate(ParentName(nInd,2))
ParentName = "NA"

deallocate(ByTmp)
deallocate(SexTmp)
deallocate(GenosR)
 
end subroutine ReadData

! ##############################################################################################

subroutine ReadParents
use Global
implicit none

integer :: i,j, dumI(2)
real :: dumR(3)
 character(len=30) :: OffName(nInd), dumC(2)

allocate(Parent(nInd,2))
Parent = 0

open(unit=103, file=trim(ParentageFileName), status="old")
read(103,*)   ! header
do i=1,nInd
    read(103,*) OffName(i), dumC(1:2), dumR(1:3), dumI(1:2), j, Parent(i,1), Parent(i,2)
    if (OffName(i) /= Id(i) .or. j/=i) then 
        call rexit("Parentage file order differs from genotype file")
    endif
enddo
 close(103)

do i=1, nInd
    if (Sex(i)==3) then
        if (ANY(Parent(:,1) == i)) then
            Sex(i) = 1
        else if (ANY(Parent(:,2) == i)) then
            Sex(i) = 2
        endif
    endif
enddo
 
end subroutine ReadParents
! ##############################################################################################

subroutine PrecalcProbs
use Global
implicit none

integer :: h,i,j,k,l,m,n
double precision :: OcA(3,3), OjA(3,3,nSnp), &
Tmp1(3), Tmp2(3,3), Tmp3(3,3,3)


!###################
! arrays of inheritance prob. conditional on 0, 1 or both parental observed genotypes
! note: all based on observed genotypes, summed over all possible actual genotypes
allocate(AHWE(3,nSnp))
allocate(OHWE(3,nSnp))
allocate(AcO(3,3,nSnp))

! Prob. observed conditional on actual
OcA(1, 1:3) = (/ 1-Er, Er/2, 0.0D0 /)   ! obs=0
OcA(2, 1:3) = (/ Er, 1-Er, Er /)    ! obs=1
OcA(3, 1:3) = (/ 0.0D0, Er/2, 1-Er /)   ! obs=2


! probabilities actual genotypes under HWE
do l=1,nSnp
    AHWE(1,l)=(1 - AF(l))**2 
    AHWE(2,l)=2*AF(l)*(1-AF(l)) 
    AHWE(3,l)=AF(l)**2 
enddo

! joined probabilities actual & observed under HWE
do l=1,nSnp
    do i=1,3    ! obs
        do j=1,3    ! act
            OjA(i, j, l) = OcA(i,j) * AHWE(j, l)
        enddo
    enddo
enddo


! marginal prob. observed genotypes
do l=1,nSnp
    do i=1,3
        OHWE(i, l) = SUM(OjA(i, 1:3, l))
    enddo
enddo


! actual genotype conditional on observed
do l=1,nSnp
    do j=1,3  !obs
        do i=1,3 !act
            AcO(i,j,l) = OjA(j,i,l)/OHWE(j,l)
        enddo
        AcO(:,j,l) = AcO(:,j,l)/SUM(AcO(:,j,l))
    enddo
enddo


! ########################
! inheritance conditional on 1 parent
allocate(AKAP(3,3,nSnp))
allocate(OKAP(3,3,nSnp))
allocate(AKOP(3,3,nSnp))
allocate(OKOP(3,3,nSnp))

do l=1,nSnp
    AKAP(1, 1:3, l) = (/ 1-AF(l), (1-AF(l))/2, 0.0D0 /)
    AKAP(2, 1:3, l) = (/ AF(l), 0.5D0, 1-AF(l) /)
    AKAP(3, 1:3, l) = (/ 0D0, AF(l)/2, AF(l) /)
enddo

do l=1,nSnp
    do i=1,3  ! obs offspring
        do j=1,3    ! act parent
            Tmp1=0
            do k=1,3    ! act offspring
                Tmp1(k) = AKAP(k,j,l) * AcO(k,i,l)
            enddo
            OKAP(i,j,l) = SUM(Tmp1)
        enddo
    enddo
enddo

do l=1,nSnp
    do i=1,3  ! act offspring
        do j=1,3    ! obs parent
            Tmp1=0
            do m=1,3    ! act parent
                Tmp1(m) = AcO(m,j,l) * AKAP(i,m,l)
            enddo
            AKOP(i,j,l) = SUM(Tmp1) 
        enddo
    enddo
enddo

do l=1,nSnp
    do i=1,3  ! obs offspring
        do j=1,3    ! obs parent
            Tmp2=0
            do k=1,3    ! act offspring
                do m=1,3    ! act parent
                    Tmp2(k,m) = AcO(m,j,l) * AKAP(k,m,l) * AcO(k,i,l)
                enddo
            enddo
            OKOP(i,j,l) = SUM(Tmp2) 
        enddo
    enddo
enddo


! #########################
! inheritance conditional on both parents

AKA2P(1,1,:) = (/ 1.0, 0.5, 0.0 /)
AKA2P(1,2,:) = (/ 0.5, 0.25, 0.0 /)
AKA2P(1,3,:) = (/ 0.0, 0.0, 0.0 /)

AKA2P(2,1,:) = (/ 0.0, 0.5, 1.0 /)
AKA2P(2,2,:) = (/ 0.5, 0.5, 0.5 /)
AKA2P(2,3,:) = (/ 1.0, 0.5, 0.0 /)

AKA2P(3,1,:) = (/ 0.0, 0.0, 0.0 /)
AKA2P(3,2,:) = (/ 0.0, 0.25, 0.5 /)
AKA2P(3,3,:) = (/ 0.0, 0.5, 1.0 /)


allocate(OKA2P(3,3,3, nSnp))
do l=1,nSnp
    do i=1,3  ! obs offspring
        do j=1,3    ! act parent 1
            do h=1,3    !act parent 2
                Tmp1=0
                do k=1,3    ! act offspring
                    Tmp1(k) = AKA2P(k,j,h) * AcO(k,i,l)
                enddo
                OKA2P(i,j,h,l) = SUM(Tmp1)
            enddo
        enddo
    enddo
enddo


allocate(OKOAP(3,3,3,nSnp))
do l=1,nSnp
    do i=1,3  ! obs offspring
        do j=1,3    ! obs parent 1
            do h=1,3    !act parent 2
                Tmp2=0
                do k=1,3    ! act offspring
                    do m=1,3    ! act parent 1
                        Tmp2(k,m) = AcO(m,j,l) * AKA2P(k,m,h) * AcO(k,i,l)
                    enddo
                enddo
                OKOAP(i,j,h,l) = SUM(Tmp2)
            enddo
        enddo
    enddo
enddo


allocate(OKO2P(3,3,3,nSnp))
do l=1,nSnp
    do i=1,3  ! obs offspring
        do j=1,3    ! obs parent 1
            do h=1,3    !obs parent 2
                Tmp3=0
                do k=1,3    ! act offspring
                    do m=1,3    ! act parent 1
                        do n=1,3    ! act parent 2
                            Tmp3(k,m,n) = AcO(m,j,l) * AcO(n,h,l) * AKA2P(k,m,n) * AcO(k,i,l)
                        enddo
                    enddo
                enddo
                OKO2P(i,j,h,l) = SUM(Tmp3)
            enddo
        enddo
    enddo
enddo

allocate(AKO2P(3,3,3,nSnp))
do l=1,nSnp
    do i=1,3  ! act offspring
        do j=1,3    ! obs parent 1
            do h=1,3    !obs parent 2
                Tmp2=0
                do m=1,3    ! act parent 1
                    do n=1,3    ! act parent 2
                        Tmp2(m,n) = AcO(m,j,l) * AcO(n,h,l) * AKA2P(i,m,n)
                    enddo
                enddo
                AKO2P(i,j,h,l) = SUM(Tmp2)
            enddo
        enddo
    enddo
enddo

allocate(AKOAP(3,3,3,nSnp))
do l=1,nSnp
    do i=1,3  ! act offspring
        do j=1,3    ! obs parent 1
            do h=1,3    !act parent 2
                Tmp1=0
                do m=1,3    ! act parent 1
                    Tmp1(m) = AcO(m,j,l) * AKA2P(i,m,h)
                enddo
                AKOAP(i,j,h,l) = SUM(Tmp1)
            enddo
        enddo
    enddo
enddo

!=================

allocate(PHS(3,3,nSnp))
allocate(PFS(3,3,nSnp))
do l=1,nSnp
    do i=1,3  ! obs offspring 1
        do j=1,3    ! obs offspring 2
            Tmp1=0
            Tmp2=0
            do m=1,3    !act shared parent 
                Tmp1(m) = OKAP(i,m,l) * OKAP(j,m,l) * AHWE(m,l)
                do h=1,3
                    Tmp2(m,h) = OKA2P(i,m,h,l) * OKA2P(j,m,h,l) * AHWE(m,l) * AHWE(h,l)
                enddo
            enddo
            PHS(i,j,l) = SUM(Tmp1) / (AHWE(i,l) * AHWE(j,l))
            PFS(i,j,l) = SUM(Tmp2) / (AHWE(i,l) * AHWE(j,l))
        enddo
    enddo
enddo

allocate(LindG(3, nSnp, nInd))  ! used when missing genotype
do l=1,nSnp
    do i=1,3
        LindG(i,l,:) = AHWE(i,l)  
    enddo
enddo

deallocate(AF)

end subroutine PrecalcProbs

! ##############################################################################################

subroutine DeAllocAll
use Global
implicit none

! allocated in ReadData
if (allocated(Sex)) deallocate(Sex)
if (allocated(BY)) deallocate(BY)
if (allocated(PairType)) deallocate(PairType)
if (allocated(nFS)) deallocate(nFS)

if (allocated(Genos)) deallocate(Genos)
if (allocated(AgeDiff)) deallocate(AgeDiff)
if (allocated(Parent)) deallocate(Parent)
if (allocated(OppHomM)) deallocate(OppHomM)
if (allocated(nS)) deallocate(nS)
if (allocated(PairID)) deallocate(PairID)
if (allocated(FSID)) deallocate(FSID)

if (allocated(SibID)) deallocate(SibID)
if (allocated(GpID)) deallocate(GpID)

if (allocated(Lind)) deallocate(Lind)
if (allocated(PairDLLR)) deallocate(PairDLLR)
if (allocated(AF)) deallocate(AF)

if (allocated(AHWE)) deallocate(AHWE)
if (allocated(OHWE)) deallocate(OHWE)
if (allocated(LLR_O)) deallocate(LLR_O)
if (allocated(LindX)) deallocate(LindX)
if (allocated(LR_parent)) deallocate(LR_parent)
if (allocated(AgePriorM)) deallocate(AgePriorM)
if (allocated(CLL)) deallocate(CLL)

if (allocated(AKAP)) deallocate(AKAP)
if (allocated(OKOP)) deallocate(OKOP)
if (allocated(AKOP)) deallocate(AKOP)
if (allocated(OKAP)) deallocate(OKAP)
if (allocated(AcO)) deallocate(AcO)
if (allocated(LR_GP)) deallocate(LR_GP)
if (allocated(LindG)) deallocate(LindG)
if (allocated(PHS)) deallocate(PHS)
if (allocated(PFS)) deallocate(PFS)

if (allocated(OKO2P)) deallocate(OKO2P)
if (allocated(OKOAP)) deallocate(OKOAP)
if (allocated(AKO2P)) deallocate(AKO2P)
if (allocated(AKOAP)) deallocate(AKOAP)
if (allocated(OKA2P)) deallocate(OKA2P)
if (allocated(DumP)) deallocate(DumP)

if (allocated(XPr)) deallocate(XPr)

if (allocated(Id)) deallocate(Id)
if (allocated(NameLH)) deallocate(NameLH)
if (allocated(ParentName)) deallocate(ParentName)

end subroutine DeAllocAll

! ##############################################################################################

! -9   NA
! 999  NA
! 888  Already assigned
! 777  impossible
! 444  not yet implemented (typically involves inbreeding)
! 222  as likely to go via opposite parent
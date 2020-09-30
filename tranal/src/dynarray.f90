!DEC$ FREEFORM
!Utility module which reshapes allocatable arrays of different dimensions and type.
module dynarray
    use precision
interface reshape  ! Generic operator for allocating and initializing or reshaping a 1D array. 
   module procedure reshapeint
   module procedure reshapeint8
   module procedure reshapereal
   module procedure reshapechar_atomname
end interface
interface reshape2  ! Generic operator for allocating and initializing or reshaping a 2D array. 
   module procedure reshapereal2
   module procedure reshapeint2
   module procedure reshapeint8_2
end interface
interface reshape3
    module procedure reshapeint3
    module procedure reshapeint8_3
    module procedure reshapereal3
end interface

contains
    subroutine reshapeint(array,newsize)
        integer(4),allocatable:: array(:),tarray(:)
        integer(4) newsize,oldsize
        if (newsize.le.0) stop 'ReshaipeInt: wrong new size!'
        if (.not. allocated(array)) then 
            allocate(array(newsize))
            array(:)=0
        elseif (size(array).lt.newsize) then
            allocate(tarray(newsize))
            oldsize=size(array)
            tarray(:)=0
            tarray(1:oldsize)=array(:)
            deallocate(array)
            allocate(array(newsize))
            array(:)=0
            array(1:oldsize)=tarray(1:oldsize)
            deallocate(tarray)                
        endif
	end subroutine reshapeint
        
        subroutine reshapeint8(array,newsize)
        integer(8),allocatable:: array(:),tarray(:)
        integer(4) newsize,oldsize
        if (newsize.le.0) stop 'ReshaipeInt: wrong new size!'
        if (.not. allocated(array)) then 
            allocate(array(newsize))
            array(:)=0
        elseif (size(array).lt.newsize) then
            allocate(tarray(newsize))
            oldsize=size(array)
            tarray(:)=0
            tarray(1:oldsize)=array(:)
            deallocate(array)
            allocate(array(newsize))
            array(:)=0
            array(1:oldsize)=tarray(1:oldsize)
            deallocate(tarray)                
        endif
	end subroutine reshapeint8

	subroutine reshapereal(array,newsize)
            real(8),allocatable:: array(:),tarray(:)
            integer(4) newsize,oldsize
            if (newsize.le.0) stop 'ReshaipeReal: New array size was not set!!!'
            if (.not. allocated(array)) then 
                allocate(array(newsize))
                array(:)=0.d0
            elseif (size(array).lt.newsize) then
                allocate(tarray(newsize))
                oldsize=size(array)
                tarray(:)=0
                tarray(1:oldsize)=array(:)
                deallocate(array)
                allocate(array(newsize))
                array(:)=0
                array(1:oldsize)=tarray(1:oldsize)
                deallocate(tarray)
            endif
	end subroutine reshapereal
        
	subroutine reshapereal2(array,newsize1,newsize2)
		real(8),allocatable:: array(:,:),tarray(:,:)
		integer(4) newsize1,newsize2
		if ((newsize1 .le. 0) .or. (newsize2 .le. 0)) stop 'ReshaipeReal2: New array size was not set!!!'
		if (.not. allocated(array)) then 
                    allocate(array(newsize1,newsize2))				
                    array(:,:)=0
		elseif ((size(array,1).lt.newsize1) .or. (size(array,2).lt.newsize2)) then				
                    allocate(tarray(max(newsize1,size(array,1)),max(newsize2,size(array,2))))
                    tarray(:,:)=0
                    tarray(1:size(array,1),1:size(array,2))=array(:,:)
                    deallocate(array)
                    allocate(array(size(tarray,1),size(tarray,2)))
                    array(:,:)=tarray(:,:)
                    deallocate(tarray)
		endif
	end subroutine reshapereal2
        
	subroutine reshapechar4(array,newsize)
		character(4),allocatable:: array(:),tarray(:)
		integer(4) newsize,oldsize
		if (newsize.le.0) stop 'ReshaipeChar4: New array size was not set!!!'
                if (.not. allocated(array)) then 
                    allocate(array(newsize))
                    array(:)=''
		elseif (size(array).lt.newsize) then
                    allocate(tarray(newsize))
                    oldsize=size(array)
                    tarray(:)=''
                    tarray(1:oldsize)=array(:)
                    deallocate(array)
                    allocate(array(newsize))
                    array(:)=tarray(:)
                    deallocate(tarray)
		endif
	end subroutine reshapechar4
        
        subroutine reshapechar_atomname(array,newsize)
		character(atom_name_len),allocatable:: array(:),tarray(:)
		integer(4) newsize,oldsize
		if (newsize.le.0) stop 'ReshaipeAtomName: New array size was not set!!!'
                if (.not. allocated(array)) then 
                    allocate(array(newsize))
                    array(:)=''
		elseif (size(array).lt.newsize) then
                    allocate(tarray(newsize))
                    oldsize=size(array)
                    tarray(:)=''
                    tarray(1:oldsize)=array(:)
                    deallocate(array)
                    allocate(array(newsize))
                    array(:)=tarray(:)
                    deallocate(tarray)
		endif
	end subroutine reshapechar_atomname
        
	subroutine reshapeint2(array,newsize1,newsize2)
		integer(4),allocatable:: array(:,:),tarray(:,:)
		integer(4) newsize1,newsize2
		if ((newsize1 .le. 0) .or. (newsize2 .le. 0)) stop 'ReshaipeInt2: New array size was not set!!!'
		if (.not. allocated(array)) then 
                    allocate(array(newsize1,newsize2))				
                    array(:,:)=0.0
		elseif ((size(array,1).lt.newsize1) .or. (size(array,2).lt.newsize2)) then                    
                    allocate(tarray(max(newsize1,size(array,1)),max(newsize2,size(array,2))))				
                    tarray(:,:)=0
                    tarray(1:size(array,1),1:size(array,2))=array(:,:)
                    deallocate(array)
                    allocate(array(size(tarray,1),size(tarray,2)))
                    array(:,:)=tarray(:,:)
                    deallocate(tarray)			
		endif
	end subroutine reshapeint2	
        
        subroutine reshapeint8_2(array,newsize1,newsize2)
		integer(8),allocatable:: array(:,:),tarray(:,:)
		integer(4) newsize1,newsize2
		if ((newsize1 .le. 0) .or. (newsize2 .le. 0)) stop 'ReshaipeInt2: New array size was not set!!!'
		if (.not. allocated(array)) then 
                    allocate(array(newsize1,newsize2))				
                    array(:,:)=0.0
		elseif ((size(array,1).lt.newsize1) .or. (size(array,2).lt.newsize2)) then                    
                    allocate(tarray(max(newsize1,size(array,1)),max(newsize2,size(array,2))))				
                    tarray(:,:)=0
                    tarray(1:size(array,1),1:size(array,2))=array(:,:)
                    deallocate(array)
                    allocate(array(size(tarray,1),size(tarray,2)))
                    array(:,:)=tarray(:,:)
                    deallocate(tarray)			
		endif
	end subroutine reshapeint8_2

        subroutine reshapeint8_3(array,newsize1,newsize2,newsize3)
		integer(8),allocatable:: array(:,:,:),tarray(:,:,:)
		integer(4) newsize1,newsize2,newsize3
		if ((newsize1 .le. 0) .or. (newsize2 .le. 0).or.(newsize3.le.0)) stop 'ReshaipeInt3: New array size was not set!!!'
		if (.not. allocated(array)) then 
                    allocate(array(newsize1,newsize2,newsize3))				
                    array=0.0
		elseif ((size(array,1).lt.newsize1) &
                    .or. (size(array,2).lt.newsize2)&
                    .or. (size(array,3).lt.newsize3)) then                    
                    allocate(tarray(max(newsize1,size(array,1)),&
                                    max(newsize2,size(array,2)),&
                                    max(newsize3,size(array,3))))				
                    tarray=0
                    tarray(1:size(array,1),1:size(array,2),1:size(array,3))=array(:,:,:)
                    deallocate(array)
                    allocate(array(size(tarray,1),size(tarray,2),size(tarray,3)))
                    array(:,:,:)=tarray(:,:,:)
                    deallocate(tarray)			
		endif
	end subroutine reshapeint8_3	
        subroutine reshapeint3(array,newsize1,newsize2,newsize3)
		integer(4),allocatable:: array(:,:,:),tarray(:,:,:)
		integer(4) newsize1,newsize2,newsize3
		if ((newsize1 .le. 0) .or. (newsize2 .le. 0).or.(newsize3.le.0)) stop 'ReshaipeInt3: New array size was not set!!!'
		if (.not. allocated(array)) then 
                    allocate(array(newsize1,newsize2,newsize3))				
                    array=0.0
		elseif ((size(array,1).lt.newsize1) &
                    .or. (size(array,2).lt.newsize2)&
                    .or. (size(array,3).lt.newsize3)) then                    
                    allocate(tarray(max(newsize1,size(array,1)),&
                                    max(newsize2,size(array,2)),&
                                    max(newsize3,size(array,3))))				
                    tarray=0
                    tarray(1:size(array,1),1:size(array,2),1:size(array,3))=array(:,:,:)
                    deallocate(array)
                    allocate(array(size(tarray,1),size(tarray,2),size(tarray,3)))
                    array(:,:,:)=tarray(:,:,:)
                    deallocate(tarray)			
		endif
	end subroutine reshapeint3	
        
        subroutine reshapereal3(array,newsize1,newsize2,newsize3)
		real(8),allocatable:: array(:,:,:),tarray(:,:,:)
		integer(4) newsize1,newsize2,newsize3
		if ((newsize1 .le. 0) .or. (newsize2 .le. 0).or.(newsize3.le.0)) stop 'ReshaipeReal3: New array size was not set!!!'
		if (.not. allocated(array)) then 
                    allocate(array(newsize1,newsize2,newsize3))				
                    array=0.0
		elseif ((size(array,1).lt.newsize1) &
                    .or. (size(array,2).lt.newsize2)&
                    .or. (size(array,3).lt.newsize3)) then                    
                    allocate(tarray(max(newsize1,size(array,1)),&
                                    max(newsize2,size(array,2)),&
                                    max(newsize3,size(array,3))))				
                    tarray=0
                    tarray(1:size(array,1),1:size(array,2),1:size(array,3))=array(:,:,:)
                    deallocate(array)
                    allocate(array(size(tarray,1),size(tarray,2),size(tarray,3)))
                    array(:,:,:)=tarray(:,:,:)
                    deallocate(tarray)			
		endif
        end subroutine reshapereal3	
        
        
end !Module

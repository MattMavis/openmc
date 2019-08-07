module utils_api_m

    use iso_c_binding

    integer(c_int), parameter, public :: MAX_ERROR_LENGTH = 1000     
    
contains

    function c_f_string(c_string, c_string_length) result(f_string)
        integer(kind=c_int), value, intent(in)                         :: c_string_length  !< The length of the c-string
        character(kind=c_char), dimension(c_string_length), intent(in) :: c_string         !< The C-String

        character(len=:), allocatable :: f_string
        integer :: i, maxlength

        maxlength = c_string_length

        do i=1, c_string_length
            if(c_string(i) == c_null_char)then
                maxlength = i - 1
                exit
            endif
        enddo

        allocate(character(len=maxlength) :: f_string)
        do i=1, maxlength
            f_string(i:i) = c_string(i)
        enddo

        f_string = trim(adjustl(f_string))    
        
    end function c_f_string
    
    function f_c_string(f_string, maxlength) result(c_string)
        character(len=*), intent(in) :: f_string    !< The fortran string
        integer(c_int), optional, intent(in) :: maxlength

        character(kind=c_char), dimension(len(f_string)+1) :: c_string
        integer :: n, i

        n = len_trim(f_string)
        if(present(maxlength))then
            n = max(0, min(len(f_string), maxlength))
        end if

        c_string(1) = c_null_char
        do i=1, n
            c_string(i) = f_string(i:i)
        enddo
        c_string(n + 1) = c_null_char

    end function f_c_string    
    
end module utils_api_m
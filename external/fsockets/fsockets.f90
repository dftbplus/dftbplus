!  F90 ISO_C_BINGING wrapper for socket communication.
!   
!  Copyright (C) 2013, Michele Ceriotti
!  Copyright (C) 2016, Bálint Aradi (changed to F2003 C-bindings & streamlined)
!   
!  Permission is hereby granted, free of charge, to any person obtaining a copy
!  of this software and associated documentation files (the "Software"), to deal
!  in the Software without restriction, including without limitation the rights
!  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!  copies of the Software, and to permit persons to whom the Software is
!  furnished to do so, subject to the following conditions:
!   
!  The above copyright notice and this permission notice shall be included in
!  all copies or substantial portions of the Software.
!   
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
!  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!  SOFTWARE.
!   
!  Contains both the functions that transmit data to the socket and read the
!  data back out again once finished, and the function which opens the socket
!  initially.
!   
!  Functions:
!
!   open_socket: Opens a socket with the required host server, socket type and
!       port number.
!
!   write_buffer: Writes a string to the socket.
!
!   read_buffer: Reads data from the socket.
!
module fsockets
  use iso_c_binding
  implicit none
  private

  public :: connect_inet_socket, connect_unix_socket, readbuffer, writebuffer, close_socket


  interface writebuffer
    module procedure writebuffer_s, writebuffer_d, writebuffer_dv, writebuffer_i
  end interface writebuffer

  interface readbuffer
    module procedure readbuffer_s
    module procedure readbuffer_dv
    module procedure readbuffer_d
    module procedure readbuffer_i
  end interface readbuffer


  ! interface bindings for c routines
  interface

    subroutine connect_inet_socket_c(sockfd, host, port)&
        & bind(C, name='fsockets_connect_inet_socket')
      import :: c_int, c_char
      integer(c_int), intent(out) :: sockfd
      character(kind=c_char), dimension(*), intent(in) :: host
      integer(c_int), value, intent(in) :: port
    end subroutine connect_inet_socket_c

    subroutine connect_unix_socket_c(sockfd, host) bind(C, name='fsockets_connect_unix_socket')
      import :: c_int, c_char
      integer(c_int), intent(out) :: sockfd
      character(kind=c_char), dimension(*), intent(in) :: host
    end subroutine connect_unix_socket_c

    subroutine writebuffer_socket_c(sockfd, data, len) bind(C, name='fsockets_writebuffer_socket')
      import :: c_int, c_ptr
      integer(c_int), value, intent(in) :: sockfd
      type(c_ptr), value, intent(in) :: data
      integer(c_int), value :: len
    end subroutine writebuffer_socket_c

    subroutine readbuffer_socket_c(sockfd, data, len) bind(C, name='fsockets_readbuffer_socket')
      import :: c_int, c_ptr
      integer(c_int), value, intent(in) :: sockfd
      type(c_ptr), value, intent(in) :: data
      integer(c_int), value, intent(in) :: len
    end subroutine readbuffer_socket_c

    subroutine close_socket_c(sockfd) bind(C, name='fsockets_close_socket')
      import :: c_int
      integer(c_int), value, intent(in) :: sockfd
    end subroutine close_socket_c

  end interface

contains

  ! Internet socket
  subroutine connect_inet_socket(sockfd, host, port)
    integer(c_int), intent(out) :: sockfd
    character(kind=c_char, len=*), intent(in) :: host
    integer(c_int), intent(in) :: port

    call connect_inet_socket_c(sockfd, trim(host) // c_null_char, port)

  end subroutine connect_inet_socket


  ! Unix socket
  subroutine connect_unix_socket(sockfd, host) 
    integer(c_int), intent(out) :: sockfd
    character(kind=c_char, len=*), intent(in) :: host

    call connect_unix_socket_c(sockfd, trim(host) // c_null_char)

  end subroutine connect_unix_socket


  ! Close the socket
  subroutine close_socket(sockfd)
    integer(c_int), intent(in) :: sockfd

    call close_socket_c(sockfd)

  end subroutine close_socket


  ! write a float (real)
  subroutine writebuffer_d(sockfd, fdata)
    integer(c_int), intent(in) :: sockfd
    real(c_double), target, intent(in) :: fdata

    call writebuffer_socket_c(sockfd, c_loc(fdata), 8_c_int)

  end subroutine writebuffer_d


  ! Write an integer 
  subroutine writebuffer_i(sockfd, fdata)
    integer(c_int), intent(in) :: sockfd
    integer(c_int), target, intent(in) :: fdata

    call writebuffer_socket_c(sockfd, c_loc(fdata), 4_c_int)

  end subroutine writebuffer_i


  ! write a string
  subroutine writebuffer_s(sockfd, fstring)
    integer(c_int), intent(in) :: sockfd
    character(kind=c_char, len=*), intent(in) :: fstring

    character(kind=c_char), target :: cstring(len(fstring))
    integer :: trimlen

    cstring(:) = transfer(fstring, cstring)
    trimlen = len_trim(fstring)
    cstring(trimlen + 1:) = ' '
    call writebuffer_socket_c(sockfd, c_loc(cstring), size(cstring, kind=c_int))
    
  end subroutine writebuffer_s


  subroutine writebuffer_dv(sockfd, fdata)
    integer(c_int), intent(in) :: sockfd
    real(c_double), target, intent(in) :: fdata(:)

    call writebuffer_socket_c(sockfd, c_loc(fdata), 8_c_int * size(fdata, kind=c_int))

  end subroutine writebuffer_dv


  subroutine readbuffer_d(sockfd, fdata)
    integer(c_int), intent(in) :: sockfd
    real(c_double), target, intent(out) :: fdata

    call readbuffer_socket_c(sockfd, c_loc(fdata), 8_c_int)

  end subroutine readbuffer_d


  subroutine readbuffer_i(sockfd, fdata)
    integer(c_int), intent(in) :: sockfd
    integer(c_int), target, intent(out) :: fdata

    call readbuffer_socket_c(sockfd, c_loc(fdata), 4_c_int)

  end subroutine readbuffer_i


  subroutine readbuffer_s(sockfd, fstring)
    integer(c_int), intent(in) :: sockfd
    character(kind=c_char, len=*), intent(out) :: fstring

    character(kind=c_char), target :: cstring(len(fstring))
    integer :: ii

    call readbuffer_socket_c(sockfd, c_loc(cstring), size(cstring, kind=c_int))
    fstring = ""
    do ii = 1, size(cstring)
      if (cstring(ii) == c_null_char) then
        exit
      end if
      fstring(ii:ii) = cstring(ii)
    end do

  end subroutine readbuffer_s


  subroutine readbuffer_dv(sockfd, fdata)
    integer(c_int), intent(in) :: sockfd
    real(c_double), target, intent(out) :: fdata(:)

    call readbuffer_socket_c(sockfd, c_loc(fdata), 8_c_int * size(fdata, kind=c_int))

  end subroutine readbuffer_dv

end module fsockets

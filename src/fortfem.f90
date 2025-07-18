module fortfem
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, fortfem!"
  end subroutine say_hello
end module fortfem

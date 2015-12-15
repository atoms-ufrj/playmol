!   This file is part of Playmol.
!
!    Playmol is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Playmol is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Playmol. If not, see <http://www.gnu.org/licenses/>.
!
!    Author: Charlles R. A. Abreu (abreu at eq.ufrj.br)
!            Applied Thermodynamics and Molecular Simulation
!            Federal University of Rio de Janeiro, Brazil
!
!   The present module is a modification of:
!
!------- -------- --------- --------- --------- --------- --------- --------- -------
! Fortran 90 function parser v1.1
!
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! * Redistribution and use in source and binary forms, with or without modification, 
! are permitted provided that the following conditions are met:
!        
! * Redistributions of source code must retain the above copyright notice, this list
! of conditions and the following disclaimer.
!        
! * Redistributions in binary form must reproduce the above copyright notice, this
! list of conditions and the following disclaimer in the documentation and/or other
! materials provided with the distribution.
!        
! * Neither the name of the copyright holder nor the names of its contributors may
! be used to endorse or promote products derived from this software without
! specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
! OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
! OF THE POSSIBILITY OF SUCH DAMAGE.
!
! The source code is available from http://fparser.sourceforge.net
!
! Roland Schmehl <roland.schmehl at alumni.uni-karlsruhe.de>
!------- -------- --------- --------- --------- --------- --------- --------- -------

module mParser

use mGlobal
use mString

implicit none

integer, parameter, private :: is = 1 ! Data type of bytecode

real(rb), parameter, private :: zero = 0.0_rb

integer(is), parameter, private :: cImmed   = 1,  &
                                   cNeg     = 2,  &
                                   cAdd     = 3,  & 
                                   cSub     = 4,  & 
                                   cMul     = 5,  & 
                                   cDiv     = 6,  & 
                                   cPow     = 7,  & 
                                   cAbs     = 8,  &
                                   cExp     = 9,  &
                                   cLog10   = 10, &
                                   cLog     = 11, &
                                   cSqrt    = 12, &
                                   cSinh    = 13, &
                                   cCosh    = 14, &
                                   cTanh    = 15, &
                                   cSin     = 16, &
                                   cCos     = 17, &
                                   cTan     = 18, &
                                   cAsin    = 19, &
                                   cAcos    = 20, &
                                   cAtan    = 21, &
                                   VarBegin = 22

character, parameter, private :: Ops(cAdd:cPow) = ["+","-", "*", "/", "^"]

character(5), parameter, private :: Funcs(cAbs:cAtan) = &
                             [character(5) :: "abs",  "exp",  "log10", "ln",  "sqrt",  &
                                              "sinh", "cosh", "tanh",  "sin", "cos",   &
                                              "tan",  "asin", "acos", "atan" ]

type tParser
  logical,              public  :: Is_Integer
  integer,              private :: NRealVar, NIntVar
  integer(is), pointer, private :: ByteCode(:)
  integer,              private :: ByteCodesize
  real(rb),    pointer, private :: Immed(:)
  logical,     pointer, private :: IsInt(:)
  integer,              private :: Immedsize
  real(rb),    pointer, private :: Stack(:)
  logical,     pointer, private :: StackInt(:)
  integer,              private :: Stacksize, StackPtr
  integer,              private :: EvalErrType
  contains
    procedure :: init => tParser_init
    procedure :: parse => tParser_parse
    procedure :: evaluate => tParser_evaluate
    procedure :: error_msg => tParser_error_msg
    procedure :: Compile => tParser_Compile
    procedure :: AddCompiledByte => tParser_AddCompiledByte
    procedure :: MathItemIndex => tParser_MathItemIndex
    procedure :: CompileSubstr => tParser_CompileSubstr
end type tParser

private :: tParser_init, tParser_parse, CheckSyntax, RemoveSpaces, Replace, &
           tParser_Compile, tParser_AddCompiledByte, tParser_CompileSubstr, LowCase

contains

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Initialize function parser
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine tParser_init( Comp )
    class(tParser), intent(inout) :: Comp
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    nullify( Comp%ByteCode, Comp%Immed, Comp%IsInt, Comp%Stack, Comp%StackInt )
  end subroutine tParser_init

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Parse function string FuncStr and compile it into bytecode
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine tParser_parse( Comp, FuncStr, RealVar, IntVar )
    class(tParser), intent(inout)        :: Comp
    character(*),   intent(in)           :: FuncStr    ! function string
    character(*),   intent(in), optional :: RealVar(:) ! Real variable names
    character(*),   intent(in), optional :: IntVar(:)  ! Integer variable names
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    character(len(FuncStr)) :: Func                    ! function string, local use
    character(sl), allocatable :: Var(:)
    integer :: Nr, Ni, ipos(len_trim(FuncStr))
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    if (present(RealVar)) then; Nr = size(RealVar); else; Nr = 0; end if
    if (present(IntVar )) then; Ni = size(IntVar);  else; Ni = 0; end if
    allocate( Var(Nr+Ni) )
    if (present(RealVar)) Var(1:Nr) = RealVar
    if (present(IntVar))  Var(Nr+1:Nr+Ni) = IntVar
    Comp%NRealVar = Nr
    Comp%NIntVar  = Ni
    nullify( Comp%ByteCode, Comp%Immed, Comp%IsInt, Comp%Stack, Comp%StackInt )
    Func = FuncStr                                  ! Local copy of function string
    call Replace("**","^ ",Func)                    ! Exponent into 1-Char. format
    call RemoveSpaces(Func,ipos)               ! Condense function string
    call CheckSyntax(Func,FuncStr,Var,ipos)
    call Comp % Compile(Func,Var)                          ! Compile into bytecode
  end subroutine tParser_parse

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Evaluate bytecode of function for the values passed in array Val(:)
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function tParser_evaluate( Comp, realval, intval ) result ( res )
    class(tParser),     intent(inout) :: Comp
    real(rb), optional, intent(in)    :: realval(:)  ! Real variable values
    integer,  optional, intent(in)    :: intval(:)   ! Integer variable values
    real(rb)                          :: res         ! result
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: IP, & ! Instruction pointer
               DP, & ! Data pointer
               SP    ! Stack pointer
    integer :: BC, Nr, Ni, ivar
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    if (present(realval)) then; Nr = size(realval); else; Nr = 0; end if
    if (present(intval))  then; Ni = size(intval);  else; Ni = 0; end if
    if (Nr /= Comp%NRealVar) &
      call error( "values of", int2str(Comp%NRealVar), "real variables are required" )
    if (Ni /= Comp%NIntVar)  &
      call error( "values of", int2str(Comp%NIntVar), "integer variables are required" )
    DP = 1
    SP = 0
    do IP = 1, Comp%ByteCodesize
      BC = Comp%ByteCode(IP)
      if ((BC > cAbs).and.(BC < VarBegin)) Comp%StackInt(SP) = .false.
      select case (BC)
        case (cImmed)
          SP = SP + 1
          Comp%Stack(SP) = Comp%Immed(DP)
          Comp%StackInt(SP) = Comp%IsInt(DP)
          DP = DP + 1
        case (cNeg) ! Negation
          Comp%Stack(SP) = -Comp%Stack(SP)
        case (cAdd) ! Addition
          Comp%StackInt(SP-1) = Comp%StackInt(SP-1) .and. Comp%StackInt(SP)
          Comp%Stack(SP-1) = Comp%Stack(SP-1) + Comp%Stack(SP)
          SP = SP - 1
        case (cSub) ! Subtraction
          Comp%StackInt(SP-1) = Comp%StackInt(SP-1) .and. Comp%StackInt(SP)
          Comp%Stack(SP-1) = Comp%Stack(SP-1) - Comp%Stack(SP)
          SP = SP - 1
        case (cMul) ! Multiplication
          Comp%StackInt(SP-1) = Comp%StackInt(SP-1) .and. Comp%StackInt(SP)
          Comp%Stack(SP-1) = Comp%Stack(SP-1) * Comp%Stack(SP)
          SP = SP - 1
        case (cDiv) ! Division
          if (Comp%Stack(SP) == zero) then
            Comp%EvalErrType = 1
            res = zero
            return
          end if
          Comp%StackInt(SP-1) = Comp%StackInt(SP-1) .and. Comp%StackInt(SP)
          if (Comp%StackInt(SP-1)) then
            Comp%Stack(SP-1) = aint(Comp%Stack(SP-1)/Comp%Stack(SP))
          else
            Comp%Stack(SP-1) = Comp%Stack(SP-1)/Comp%Stack(SP)
          end if
          SP = SP - 1
        case (cPow) ! Powering
          Comp%StackInt(SP-1) = Comp%StackInt(SP-1) .and. Comp%StackInt(SP)
          if (Comp%StackInt(SP-1)) then
            Comp%Stack(SP-1) = int(Comp%Stack(SP-1))**int(Comp%Stack(SP))
          else if (Comp%StackInt(SP)) then
            Comp%Stack(SP-1) = Comp%Stack(SP-1)**int(Comp%Stack(SP))
          else
            Comp%Stack(SP-1) = Comp%Stack(SP-1)**Comp%Stack(SP)
          end if
          SP = SP - 1
        case (cAbs) ! Absolute value
          Comp%Stack(SP) = abs(Comp%Stack(SP))
        case (cExp) ! Exponentiation
          Comp%Stack(SP) = exp(Comp%Stack(SP))
        case (cLog10) ! Decimal logarithm
          if (Comp%Stack(SP) <= zero) then
            Comp%EvalErrType = 3
            res = zero
            return
          end if
          Comp%Stack(SP) = log10(Comp%Stack(SP))
        case (cLog) ! Natural logarithm
          if (Comp%Stack(SP) <= zero) then
            Comp%EvalErrType = 3
            res = zero
            return
          end if
          Comp%Stack(SP) = log(Comp%Stack(SP))
        case (cSqrt) ! Square root
          if (Comp%Stack(SP) < zero) then
            Comp%EvalErrType = 3
            res = zero
            return
          end if
          Comp%Stack(SP) = sqrt(Comp%Stack(SP))
        case (cSinh) ! Hyperbolic sine
          Comp%Stack(SP) = sinh(Comp%Stack(SP))
        case (cCosh) ! Hyperbolic cosine
          Comp%Stack(SP) = cosh(Comp%Stack(SP))
        case (cTanh) ! Hyberbolic tangent
          Comp%Stack(SP) = tanh(Comp%Stack(SP))
        case (cSin)  ! Sine
          Comp%Stack(SP) = sin(Comp%Stack(SP))
        case (cCos)  ! Cosine
          Comp%Stack(SP) = cos(Comp%Stack(SP))
        case (cTan)  ! Tangent
          Comp%Stack(SP) = tan(Comp%Stack(SP))
        case (cAsin) ! Arcsine
          if ((Comp%Stack(SP) < -1.0_rb).or.(Comp%Stack(SP) > 1.0_rb)) then
            Comp%EvalErrType = 4
            res = zero
            return
          end if
          Comp%Stack(SP) = asin(Comp%Stack(SP))
        case (cAcos) ! Arccosine
          if ((Comp%Stack(SP) < -1.0_rb).or.(Comp%Stack(SP) > 1.0_rb)) then
            Comp%EvalErrType = 4
            res = zero
            return
          end if
          Comp%Stack(SP) = acos(Comp%Stack(SP))
        case (cAtan) ! Arctangent
          Comp%Stack(SP) = atan(Comp%Stack(SP))
        case default
          SP = SP + 1
          ivar = BC - VarBegin + 1
          Comp%StackInt(SP) = ivar > Nr
          if (Comp%StackInt(SP)) then
            Comp%Stack(SP) = intval(ivar-Nr)
          else
            Comp%Stack(SP) = realval(ivar)
          end if
        end select
     end do
     Comp%EvalErrType = 0
     res = Comp%Stack(1)
     Comp%Is_Integer = Comp%StackInt(1)
  end function tParser_evaluate

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return error message
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function tParser_error_msg( Comp ) result (msg)
    class(tParser), intent(in) :: Comp
    character(*), parameter :: m(4) = [ "Division by zero                ", &
                                        "Argument of SQRT negative       ", &
                                        "Argument of LOG negative        ", &
                                        "Argument of ASIN or ACOS illegal" ]
    character(len(m))       :: msg
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    if ((Comp%EvalErrType < 1) .or. (Comp%EvalErrType > size(m))) then
       msg = ""
    else
       msg = m(Comp%EvalErrType)
    end if
  end function tParser_error_msg

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Check syntax of function string,  returns 0 if syntax is ok
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine CheckSyntax( Func, FuncStr, Var, ipos )
    character(*), intent(in) :: Func         ! function string without spaces
    character(*), intent(in) :: FuncStr      ! Original function string
    character(*), intent(in) :: Var(:)       ! Array with variable names
    integer,      intent(in) :: ipos(:)
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer(is) :: n
    character   :: c
    real(rb)    :: r
    logical     :: err
    integer     :: ParCnt, & ! Parenthesis counter
                   j,ib,in,lFunc
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    j = 1
    ParCnt = 0
    lFunc = len_trim(Func)
    step: do
      if (j > lFunc) call ParseErrMsg(j, FuncStr)
      c = Func(j:j)
      ! Check for valid operand (must appear)
      if ((c == "-") .or. (c == "+")) then           ! Check for leading - or +
        j = j + 1
        if (j > lFunc) &
          call ParseErrMsg(j, FuncStr, "Missing operand")
        c = Func(j:j)
        if (any(c == Ops)) &
          call ParseErrMsg(j, FuncStr, "Multiple operators")
      end if
      n = MathFunctionIndex(Func(j:))
      if (n > 0) then                                ! Check for math function
        j = j + len_trim(Funcs(n))
        if (j > lFunc) &
          call ParseErrMsg(j, FuncStr, "Missing function argument")
        c = Func(j:j)
        if (c /= "(") &
          call ParseErrMsg(j, FuncStr, "Missing opening parenthesis")
      end if
      if (c == "(") then                             ! Check for opening parenthesis
        ParCnt = ParCnt + 1
        j = j + 1
        cycle step
      end if
      if (scan(c,"0123456789.") > 0) then            ! Check for number
        r = realNum(Func(j:),ib,in,err)
        if (err) call ParseErrMsg(j, FuncStr, &
                                  "Invalid number format: "//Func(j+ib-1:j+in-2))
        j = j + in - 1
        if (j > lFunc) exit
        c = Func(j:j)
      else                                           ! Check for variable
        n = VariableIndex (Func(j:),Var,ib,in)
        if (n == 0) call ParseErrMsg(j, FuncStr, &
                                      "Invalid element: "//Func(j+ib-1:j+in-2))
        j = j + in - 1
        if (j > lFunc) exit
        c = Func(j:j)
      end if
      do while (c == ")")                            ! Check for closing parenthesis
        ParCnt = ParCnt - 1
        if (ParCnt < 0) &
          call ParseErrMsg(j, FuncStr, "Mismatched parenthesis")
        if (Func(j-1:j-1) == "(") &
          call ParseErrMsg(j-1, FuncStr, "Empty parentheses")
        j = j + 1
        if (j > lFunc) exit
        c = Func(j:j)
      end do
      ! Now, we have a legal operand: A legal operator or end of string must follow
      if (j > lFunc) exit
      if (any(c == Ops)) then                        ! Check for multiple operators
        if (j+1 > lFunc) &
          call ParseErrMsg(j, FuncStr)
        if (any(Func(j+1:j+1) == Ops)) &
          call ParseErrMsg(j+1, FuncStr, "Multiple operators")
      else                                           ! Check for next operand
        call ParseErrMsg(j, FuncStr, "Missing operator")
      end if
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      ! Now, we have an operand and an operator: the next loop will check for another 
      ! operand (must appear)
      !-- -------- --------- --------- --------- --------- --------- --------- -------
      j = j + 1
    end do step
    if (ParCnt > 0) call ParseErrMsg(j, FuncStr, "Missing )")

    contains

      ! Print error message and terminate program
      subroutine ParseErrMsg( j, FuncStr, Msg )
        integer,                intent(in) :: j
        character(*),           intent(in) :: FuncStr       ! Original function string
        character(*), optional, intent(in) :: Msg
        !-------- --------- --------- --------- --------- --------- --------- -------
        call writeln( "Error: wrong syntax of function string" )
        if (present(Msg)) then
           call writeln( ">", Msg )
        end if
        call writeln( ">", FuncStr )
        call writeln( ">", repeat(" ",ipos(j)-1)//"^" )
        stop
      end subroutine ParseErrMsg

  end subroutine CheckSyntax

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return operator index
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function OperatorIndex( c ) result (n)
    character(1), intent(in) :: c
    integer(is)              :: n
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer (is) :: j
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    do j = cAdd, cPow
       if (c == Ops(j)) then
          n = j
          exit
       end if
    end do
  end function OperatorIndex

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return index of math function beginnig at 1st position of string str
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function MathFunctionIndex( str ) result (n)
    character(*), intent(in) :: str
    integer(is)              :: n
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer(is) :: j
    integer     :: k
    character(len(Funcs)) :: fun
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    do j = cAbs, cAtan                                   ! Check all math functions
       k = min(len_trim(Funcs(j)), len(str))   
       call LowCase (str(1:k), fun)
       if (fun == Funcs(j)) then                         ! Compare lower case letters
          n = j                                          ! Found a matching function
          exit
       end if
    end do
  end function MathFunctionIndex

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return index of variable at begin of string str (returns 0 if no variable found)
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function VariableIndex( str, Var, ibegin, inext ) result (n)
    character(*),      intent(in)  :: str       ! String
    character(*),      intent(in)  :: Var(:)    ! Array with variable names
    integer, optional, intent(out) :: ibegin, & ! Start position of variable name
                                      inext     ! Position of character after name
    integer(is)                    :: n         ! index of variable
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: j,ib,in,lstr
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    lstr = len_trim(str)
    if (lstr > 0) then
       do ib = 1, lstr                          ! Search for first character in str
          if (str(ib:ib) /= " ") exit           ! When lstr>0 at least 1 char in str
       end do                        
       do in = ib, lstr                         ! Search for name terminators
          if (scan(str(in:in),"+-*/^) ") > 0) exit
       end do
       do j = 1, size(Var)
          if (str(ib:in-1) == Var(j)) then                     
             n = j                              ! Variable name found
             exit
          end if
       end do
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
  end function VariableIndex

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Remove Spaces from string, remember positions of characters in old string
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine RemoveSpaces( str, ipos )
    character(*), intent(inout) :: str
    integer,      intent(inout) :: ipos(:)
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: k, lstr
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    lstr = len_trim(str)
    ipos = [(k,k=1,lstr)]
    k = 1
    do while (str(k:lstr) /= " ")                             
       if (str(k:k) == " ") then
          str(k:lstr)  = str(k+1:lstr)//" "    ! Move 1 character to left
          ipos(k:lstr) = [ipos(k+1:lstr), 0]   ! Move 1 element to left
          k = k-1
       end if
       k = k+1
    end do
  end subroutine RemoveSpaces

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Replace ALL appearances of character set ca in string str by character set cb
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine Replace( ca, cb, str )
    character(*),       intent(in) :: ca
    character(len(ca)), intent(in) :: cb                ! len(ca) must be len(cb)
    character(*),    intent(inout) :: str
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer  :: j, lca
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    lca = len(ca)
    do j = 1, len_trim(str)-lca+1
       if (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
    end do
  end subroutine Replace

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Compile i-th function string F into bytecode
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine tParser_Compile( Comp, F, Var )
    class(tParser), intent(inout) :: Comp
    character(*),   intent(in)    :: F            ! function string
    character(*),   intent(in)    :: Var(:)       ! Array with variable names
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: istat
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    if (associated(Comp%ByteCode)) deallocate ( Comp%ByteCode, &
                                                Comp%Immed,    &
                                                Comp%Stack     )
    Comp%ByteCodesize = 0
    Comp%Immedsize    = 0
    Comp%Stacksize    = 0
    Comp%StackPtr     = 0
    call Comp % CompileSubstr(F,1,len_trim(F),Var)     ! Compile string to determine size
    allocate( Comp%ByteCode(Comp%ByteCodesize), &
              Comp%Immed(Comp%Immedsize),       &
              Comp%IsInt(Comp%Immedsize),       &
              Comp%Stack(Comp%Stacksize),       &
              Comp%StackInt(Comp%Stacksize),    &
              STAT = istat                      )
    if (istat /= 0) then
      call error( "Memmory allocation for byte code failed" )
    else
      Comp%ByteCodesize = 0
      Comp%Immedsize    = 0
      Comp%Stacksize    = 0
      Comp%StackPtr     = 0
      call Comp % CompileSubstr(F,1,len_trim(F),Var)  ! Compile string into bytecode
    end if
  end subroutine tParser_Compile

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Add compiled byte to bytecode
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine tParser_AddCompiledByte( Comp, b )
    class(tParser), intent(inout) :: Comp
    integer(is),    intent(in)    :: b                   ! Value of byte to be added
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    Comp%ByteCodesize = Comp%ByteCodesize + 1
    if (associated(Comp%ByteCode)) Comp%ByteCode(Comp%ByteCodesize) = b
  end subroutine tParser_AddCompiledByte

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! return math item index, if item is real number, enter it into Comp-structure
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function tParser_MathItemIndex( Comp, F, Var ) result( n )
    class(tParser), intent(inout) :: Comp
    character(*),   intent(in)    :: F         ! function substring
    character(*),   intent(in)    :: Var(:)    ! Array with variable names
    integer(is)                   :: n         ! Byte value of math item
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    if (scan(F(1:1),"0123456789.") > 0) then            ! Check for begin of a number
      Comp%Immedsize = Comp%Immedsize + 1
      if (associated(Comp%Immed)) then
        Comp%Immed(Comp%Immedsize) = RealNum(F)
        Comp%IsInt(Comp%Immedsize) = IsIntNum(F)
      end if
      n = cImmed
    else                                                ! Check for a variable
      n = VariableIndex(F, Var)
      if (n > 0) n = VarBegin + n - 1
    end if
  end function tParser_MathItemIndex

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function CompletelyEnclosed( F, b, e ) result (res)
    character(*), intent(in) :: F                   ! function substring
    integer,      intent(in) :: b,e                 ! First and last pos. of substring
    logical                  :: res
    !---- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: j, k
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    res = .false.
    if (F(b:b) == "(" .and. F(e:e) == ")") then
      k = 0
      do j = b+1, e-1
        if (F(j:j) == "(") then
          k = k+1
        else if (F(j:j) == ")") then
          k = k-1
        end if
        if (k < 0) exit
      end do
      if (k == 0) res = .true.                      ! All opened parenthesis closed
    end if
  end function CompletelyEnclosed

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Compile function string F into bytecode
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  recursive subroutine tParser_CompileSubstr( Comp, F, b, e, Var )
    class(tParser), intent(inout) :: Comp
    character(*),   intent(in)    :: F         ! function substring
    integer,        intent(in)    :: b,e       ! Begin and end position substring
    character(*),   intent(in)    :: Var(:)    ! Array with variable names
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer(is) :: n
    integer     :: b2,j,k,io
    character(*), parameter :: calpha = "abcdefghijklmnopqrstuvwxyz"// &
                                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    ! Check for special cases of substring
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    if (F(b:b) == "+") then                              ! case 1: F(b:e) = "+..."
!     write(*,*)"1. F(b:e) = "+...""
      call Comp % CompileSubstr( F, b+1, e, Var )
      return
    else if (CompletelyEnclosed(F, b, e)) then           ! case 2: F(b:e) = "(...)"
!     write(*,*)"2. F(b:e) = "(...)""
      call Comp % CompileSubstr( F, b+1, e-1, Var )
      return
    else if (scan(F(b:b),calpha) > 0) then        
      n = MathFunctionIndex (F(b:e))
      if (n > 0) then
        b2 = b + index(F(b:e),"(") - 1
        if (CompletelyEnclosed(F, b2, e)) then         ! case 3: F(b:e) = "fcn(...)"
!         write(*,*)"3. F(b:e) = "fcn(...)""
          call Comp % CompileSubstr( F, b2+1, e-1, Var )
          call Comp % AddCompiledByte( n )
          return
        end if
      end if
    else if (F(b:b) == "-") then
      if (CompletelyEnclosed(F, b+1, e)) then           ! case 4: F(b:e) = "-(...)"
!       write(*,*)"4. F(b:e) = "-(...)""
        call Comp % CompileSubstr( F, b+2, e-1, Var )
        call Comp % AddCompiledByte( cNeg )
        return
      else if (scan(F(b+1:b+1),calpha) > 0) then
        n = MathFunctionIndex(F(b+1:e))
        if (n > 0) then
          b2 = b + index(F(b+1:e),"(")
          if (CompletelyEnclosed(F, b2, e)) then      ! case 5: F(b:e) = "-fcn(...)"
!           write(*,*)"5. F(b:e) = "-fcn(...)""
            call Comp % CompileSubstr( F, b2+1, e-1, Var )
            call Comp % AddCompiledByte( n )
            call Comp % AddCompiledByte( cNeg )
            return
          end if
        end if
      end if
    end if

    ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
    do io = cAdd, cPow                                       ! Increasing priority +-*/^
      k = 0
      do j = e, b, -1
        if (F(j:j) == ")") then
           k = k + 1
        else if (F(j:j) == "(") then
           k = k - 1
        end if
        if ((k == 0) .and. (F(j:j) == Ops(io)) .and. (IsBinaryOp(j, F))) then
          if (any(F(j:j) == Ops(cMul:cPow)) .and. (F(b:b) == "-")) then ! case 6: F(b:e) = "-...Op..." with Op > -
!           write(*,*)"6. F(b:e) = "-...Op..." with Op > -"
            call Comp % CompileSubstr( F, b+1, e, Var )
            call Comp % AddCompiledByte( cNeg)
            return                 
          else                                                        ! case 7: F(b:e) = "...BinOp..."
!           write(*,*)"7. Binary operator",F(j:j)
            call Comp % CompileSubstr( F, b, j-1, Var )
            call Comp % CompileSubstr( F, j+1, e, Var )
            call Comp % AddCompiledByte( OperatorIndex(Ops(io)) )
            Comp%StackPtr = Comp%StackPtr - 1
            return
          end if
        end if
      end do
    end do
    ! Check for remaining items, i.e. variables or explicit numbers
    b2 = b
    if (F(b:b) == "-") b2 = b2 + 1
    n = Comp % MathItemIndex( F(b2:e), Var )
!   write(*,*)"8. AddCompiledByte ",n
    call Comp % AddCompiledByte( n)
    Comp%StackPtr = Comp%StackPtr + 1
    if (Comp%StackPtr > Comp%Stacksize) Comp%Stacksize = Comp%Stacksize + 1
    if (b2 > b) call Comp % AddCompiledByte( cNeg )
  end subroutine tParser_CompileSubstr

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Check if operator F(j:j) in string F is binary operator
  ! Special cases already covered elsewhere:              (that is corrected in v1.1)
  ! - operator character F(j:j) is first character of string (j=1)
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function IsBinaryOp( j, F ) result (res)
    integer,      intent(in) :: j                       ! Position of Operator
    character(*), intent(in) :: F                       ! String
    logical                  :: res                     ! result
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: k
    logical :: Dflag,Pflag
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    res=.true.
    if ((F(j:j) == "+") .or. (F(j:j) == "-")) then       ! Plus or minus sign:
      if (j == 1) then                                   ! - leading unary operator ?
        res = .false.
      else if (scan(F(j-1:j-1),"+-*/^(") > 0) then       ! - other unary operator ?
        res = .false.
      else if (scan(F(j+1:j+1),"0123456789") > 0 .and. & ! - in exponent of real number ?
               scan(F(j-1:j-1),"eEdD")       > 0) then
        Dflag=.false.; Pflag=.false.
        k = j - 1
        do while (k > 1)                                 !   step to the left in mantissa 
          k = k - 1
          if (scan(F(k:k),"0123456789") > 0) then
            Dflag=.true.
          else if (F(k:k) == ".") then
            if (Pflag) then
              exit                                       !   * exit: 2nd appearance of "."
            else
              Pflag=.true.                               !   * mark 1st appearance of "."
            end if
          else
            exit                                         !   * all other characters
          end if
        end do
        if (Dflag .and. (k == 1 .or. scan(F(k:k),"+-*/^(") > 0)) res = .false.
      end if
    end if
  end function IsBinaryOp

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function realNum( str, ibegin, inext, error ) result (res)
    character(*),      intent(in)  :: str           ! String
    integer, optional, intent(out) :: ibegin,     & ! Start position of real number
                                      inext         ! 1st character after real number
    logical, optional, intent(out) :: error         ! Error flag
    real(rb)                       :: res           ! real number
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer  :: ib,in,istat
    logical  :: Bflag,               & ! .T. at begin of number in str
                InMan,               & ! .T. in mantissa of number
                Pflag,               & ! .T. after 1st "." encountered
                Eflag,               & ! .T. at exponent identifier "eEdD"
                InExp,               & ! .T. in exponent of number
                DInMan,              & ! .T. if at least 1 digit in mant.
                DInExp,              & ! .T. if at least 1 digit in exp.
                err                    ! Local error flag
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    Bflag = .true.
    InMan = .false.
    Pflag = .false.
    Eflag = .false.
    InExp = .false.
    DInMan = .false.
    DInExp = .false.
    ib = 1
    in = 1
    do while (in <= len_trim(str))
      select case (str(in:in))
      case (" ")                               ! Only leading blanks permitted
        ib = ib + 1
        if (InMan .or. Eflag .or. InExp) exit
      case ("+","-")                           ! Permitted only
        if (Bflag) then           
          InMan = .true.
          Bflag = .false.                      ! - at beginning of mantissa
        else if (Eflag) then               
          InExp = .true.
          Eflag = .false.                      ! - at beginning of exponent
        else
          exit                                 ! - otherwise stop
        end if
      case ("0":"9")                           ! Mark
        if (Bflag) then           
          InMan = .true.
          Bflag = .false.                      ! - beginning of mantissa
        else if (Eflag) then               
           InExp = .true.
           Eflag = .false.                     ! - beginning of exponent
        end if
        if (InMan) DInMan=.true.               ! Mantissa contains digit
        if (InExp) DInExp=.true.               ! Exponent contains digit
      case (".")
        if (Bflag) then
          Pflag = .true.                       ! - mark 1st appearance of "."
          InMan = .true.
          Bflag = .false.                      !   mark beginning of mantissa
        else if (InMan .and..not.Pflag) then
          Pflag = .true.                       ! - mark 1st appearance of "."
        else
           exit                                ! - otherwise stop
        end if
      case ("e","E","d","D")                   ! Permitted only
        if (InMan) then
           Eflag = .true.
           InMan = .false.                     ! - following mantissa
        else
           exit                                ! - otherwise stop
        end if
      case default
        exit                                   ! stop at all other characters
      end select
      in = in + 1
    end do
    err = (ib > in-1) .or. (.not.DInMan) .or. ((Eflag.or.InExp).and..not.DInExp)
    if (err) then
      res = zero
    else
      read(str(ib:in-1),*,iostat=istat) res
      err = istat /= 0
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
    if (present(error))  error  = err
  end function realNum

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Check integer number format: [blanks][+|-][nnn][blanks]
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  function IsIntNum( str ) result( res )
    character(*),      intent(in)  :: str           ! String
    logical                        :: res           ! True if str is an integer
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    character(len_trim(str)) :: C
    integer :: k, n
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    C = adjustl(str)
    n = len_trim(C)
    if (n == 0) then
      res = .false.
    else if (n == 1) then
      res = scan(C(1:1),"0123456789") /= 0
    else
      res = scan(C(1:1),"+-0123456789") /= 0
      k = 1
      do while (res .and. (k < n))
        k = k + 1
        res = scan(C(k:k),"+-0123456789") /= 0
      end do
    end if
  end function IsIntNum

  !----- -------- --------- --------- --------- --------- --------- --------- -------
  ! Transform upper case letters in str1 into lower case letters, result is str2
  !----- -------- --------- --------- --------- --------- --------- --------- -------
  subroutine LowCase( str1, str2 )
    character(*), intent(in)  :: str1
    character(*), intent(out) :: str2
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    integer :: j, k
    character(*), parameter :: lc = "abcdefghijklmnopqrstuvwxyz"
    character(*), parameter :: uc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    !--- -------- --------- --------- --------- --------- --------- --------- -------
    str2 = str1
    do j = 1, len_trim(str1)
      k = index(uc,str1(j:j))
      if (k > 0) str2(j:j) = lc(k:k)
    end do
  end subroutine LowCase

end module mParser

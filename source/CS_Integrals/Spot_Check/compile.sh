EXE="x_test_int.exe" 

#OPT="-O0 -g -Mprof=func -Mpreprocess -C "
OPT="-O0 -C -traceback -g "
LAPACK="-llapack"
BLAS="-lblas"
FC77="ifort "

rm  $EXE *.o

echo "COMPILING...."

set echo

$FC77 $OPT  ../../AC_Integrals/*.f ../*.f *.f  -o  $EXE

unset echo

echo "DONE!"


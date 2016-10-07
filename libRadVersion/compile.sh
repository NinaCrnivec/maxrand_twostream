rm -f run_twomaxrnd *.o

DL=$HOME/libRadtran/

make -j -C $DL libRadtran_c

gcc -c main_twomaxrnd.c -DF77_FUNC\(name,NAME\)=name\ ##\ _ -I$DL/libsrc_c -I$DL/libsrc_f

gfortran -o run_twomaxrnd main_twomaxrnd.o -L$DL/lib/ -lRadtran_c -lRadtran_f



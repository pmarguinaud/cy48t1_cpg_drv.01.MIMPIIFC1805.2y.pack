#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=arp
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --exclusive
#SBATCH -p ndl
#SBATCH --verbose
#SBATCH --no-requeue

set -x

# Environment variables

ulimit -s unlimited
export OMP_STACKSIZE=4G
export KMP_STACKSIZE=4G
export KMP_MONITOR_STACKSIZE=1G
export DR_HOOK=1
export DR_HOOK_IGNORE_SIGNALS=-1
export DR_HOOK_OPT=prof
export EC_PROFILE_HEAP=0
export EC_MPI_ATEXIT=0

export PATH=$HOME/bin:$HOME/benchmf1709/scripts:$PATH

# Directory where input data is stored

DATA=/scratch/work/marguina/benchmf1709-data

# Change to a temporary directory

export workdir=/scratch/work/marguina

if [ "x$SLURM_JOBID" != "x" ]
then
export TMPDIR=$workdir/tmp/arp.$SLURM_JOBID
else
export TMPDIR=$workdir/tmp/arp.$$
fi

mkdir -p $TMPDIR

cd $TMPDIR


#for K in 0 1
for K in 1 4
do

mkdir -p $K
cd $K


# Choose your test case resolution

#GRID=t1798
#GRID=t0798
 GRID=t0031

# Choose a pack

 PACK=/home/gmap/mrpm/marguina/pack/cy48t1_cpg_drv.01.MIMPIIFC1805.2y.pack

# Copy data to $TMPDIR

for f in $DATA/arp/grid/$GRID/* $DATA/arp/code/43oops/data/* 
do
  if [ -f "$f" ]
  then
    \rm -f $(basename $f)
    ln $f .
  fi
done

for f in $PACK/data/*
do
  b=$(basename $f)
  \rm -f $b
  ln -s $f $b
done

for f in $DATA/arp/code/43oops/naml/*
do
  cp $f .
  chmod 644 $(basename $f)
done

# Set the number of nodes, tasks, threads for the model

NNODE_FC=1
NTASK_FC=1
NOPMP_FC=1

# Set the number of nodes, tasks, threads for the IO server

NNODE_IO=0
NTASK_IO=32
NOPMP_IO=4

let "NPROC_FC=$NNODE_FC*$NTASK_FC"
let "NPROC_IO=$NNODE_IO*$NTASK_IO"

# Set forecast term; reduce it for debugging

STOP=6

# Modify namelist

xpnam --delta="
&NAMRIP
! CSTOP='h$STOP',
  CSTOP='t10',
  TSTEP=240,
/
&NAMARG
  CUSTOP=-,
  UTSTEP=-,
/
&NAMCVER
  NVSCH=-,
  NVFE_TYPE=1,
/
&NAMTRAJ
/
&NAMCT0
  LGRIB_API=.FALSE.,
/
&NAETLDIAG
/
&NAMDVISI
/
&NAMSATSIM
/
&NAMPHY0
  GREDDRS=-,
/
&NAMNORGWD
/
&NAMMETHOX
/
&NAMNUDGLH
/
&NAMSPP
/
&NAMFA
  NVGRIB=0,
/
&NAMDPRECIPS
/
" --inplace fort.4

xpnam --delta="
&NAMIO_SERV
  NPROC_IO=$NPROC_IO,
/
&NAMPAR0
  NPROC=$NPROC_FC,
  $(square $NPROC_FC)
/
&NAMPAR1
  NSTRIN=$NPROC_FC,
/
" --inplace fort.4

# Set up grib_api environment

grib_api_prefix=$(ldd $PACK/bin/MASTERODB  | perl -ne ' 
   next unless (m/libeccodes_f90/o); 
   my ($path) = (m/=>\s+(\S+)/o); 
   use File::Basename; 
   print &dirname (&dirname ($path)), "\n" ')
export GRIB_DEFINITION_PATH=$PWD/extra_grib_defs:$grib_api_prefix/share/definitions:$grib_api_prefix/share/eccodes/definitions
export GRIB_SAMPLES_PATH=$grib_api_prefix/ifs_samples/grib1:$grib_api_prefix/share/eccodes/ifs_samples/grib1

# Change NPROMA

if [ 1 -eq 1 ]
then
xpnam --delta="
&NAMDIM
  NPROMA=-32,
/
" --inplace fort.4
fi

xpnam --delta="
$(cat $PACK/nam/APLPAR_NEW.nam)
" --inplace fort.4


xpnam --delta="
&NAMARPHY
  LAPL_ARPEGE=.TRUE.
/
" --inplace fort.4


ls -lrt

cat fort.4

# Run the model; use your mpirun

if [ "x$K" = "x0" ]
then
  export INPART=0
  export PARALLEL=0
  export PERSISTENT=0
  export OPENACC=0
elif [ "x$K" = "x1" ]
then
  export INPART=1
  export PARALLEL=1
  export PERSISTENT=1
  export OPENACC=0
elif [ "x$K" = "x2" ]
then
  export INPART=0
  export PARALLEL=0
  export PERSISTENT=1
  export OPENACC=0
elif [ "x$K" = "x3" ]
then
  export INPART=1
  export PARALLEL=0
  export PERSISTENT=1
  export OPENACC=0
elif [ "x$K" = "x4" ]
then
  export INPART=1
  export PARALLEL=1
  export PERSISTENT=1
  export OPENACC=1
else
  unset INPART
  unset PARALLEL
  unset PERSISTENT
  unset OPENACC
fi


pack=$PACK
pack=/home/gmap/mrpm/marguina/pack/48t1_cpg_drv.01.PGI217.cpu0
BIN=$pack/bin/MASTERODB

export DEBUG=0


if [ "x$DEBUG" = "x1" ]
then

(
/opt/softs/mpiauto/mpiauto --verbose --wrap --wrap-stdeo --nouse-slurm-mpi --prefix-mpirun '/usr/bin/time -f "time=%e"' \
    --nnp $NTASK_FC --nn $NNODE_FC --openmp $NOPMP_FC -- $BIN \
 -- --nnp $NTASK_IO --nn $NNODE_IO --openmp $NOPMP_IO -- $BIN 
) &


for i in 0 
do
while [ True ]
do
  if [ -f "pid.$i.txt" ]
  then
    break
  fi  
  sleep 1
done
done
set -e


gdb -ex 'shell rm pid.0.txt' $BIN $(cat pid.0.txt) 

rm pid.*.txt

wait

else

export UTIL="COPY=MODEL,GEOMETRY"

/opt/softs/mpiauto/mpiauto --verbose --wrap --wrap-stdeo --nouse-slurm-mpi --prefix-mpirun '/usr/bin/time -f "time=%e"' \
    --nnp $NTASK_FC --nn $NNODE_FC --openmp $NOPMP_FC -- $BIN \
 -- --nnp $NTASK_IO --nn $NNODE_IO --openmp $NOPMP_IO -- $BIN 


fi



ls -lrt

cd ..

done


#diffNODE.001_01 0/NODE.001_01 2/NODE.001_01
#diffNODE.001_01 0/NODE.001_01 1/NODE.001_01
diffNODE.001_01 1/NODE.001_01 4/NODE.001_01



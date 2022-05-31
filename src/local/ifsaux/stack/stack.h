#ifndef _STACK_N
#define _STACK_N

#define USE_STACK

#ifdef USE_STACK

#define temp(t, n, s) t, DIMENSION s :: n; POINTER (IP_##n##_, n)

USE ABOR1_ACC_MOD

#define alloc(n) IP_##n##_=YLSTACK%L;YLSTACK%L=YLSTACK%L+JPRB*SIZE(n);IF(YLSTACK%L>YLSTACK%U)CALL ABOR1_ACC(__FILE__)

#else

#define temp(t, n, s) t, DIMENSION s :: n
#define alloc(n)

#endif

#endif

#ifdef ESMC_RCS_HEADER
"$Id: conf.h,v 1.3 2002/12/16 22:43:37 nscollins Exp $"
"Defines the configuration for this machine"
#endif

#if 0
"i have NOT checked any of these values.  i just copied"
" them from the linux directory"
#endif

#if !defined(INCLUDED_CONF_H)
#define INCLUDED_CONF_H

#define PARCH_mac_osx
#define ESMF_ARCH_NAME "mac_osx"

#define ESMC_HAVE_LIMITS_H
#define ESMC_HAVE_PWD_H 
#define ESMC_HAVE_MALLOC_H 
#define ESMC_HAVE_STRING_H 
#define ESMC_HAVE_GETDOMAINNAME
#define ESMC_HAVE_DRAND48 
#define ESMC_HAVE_UNAME 
#define ESMC_HAVE_UNISTD_H 
#define ESMC_HAVE_SYS_TIME_H 
#define ESMC_HAVE_STDLIB_H

#define ESMC_SUBSTITUTE_CTRL_CHARS 1

#define ESMC_POINTER_SIZE 4
#undef ESMC_HAVE_OMP_THREADS 

#define ESMC_HAVE_MPI 1

#define ESMC_HAVE_READLINK
#define ESMC_HAVE_MEMMOVE
#define ESMC_HAVE_TEMPLATED_COMPLEX

#define ESMC_HAVE_DOUBLE_ALIGN_MALLOC
#define ESMC_HAVE_MEMALIGN
#define ESMC_HAVE_SYS_RESOURCE_H
#define ESMC_SIZEOF_VOIDP 4
#define ESMC_SIZEOF_INT 4
#define ESMC_SIZEOF_DOUBLE 8

#if defined(fixedsobug)
#define ESMC_USE_DYNAMIC_LIBRARIES 1
#define ESMC_HAVE_RTLD_GLOBAL 1
#endif

#define ESMC_HAVE_SYS_UTSNAME_H

#define ESMF_IS_32BIT_MACHINE 1

#endif

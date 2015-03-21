//this is just a copy of lapacke_mangling_with_flags.h, not really sure why I needed to do that
#ifndef LAPACK_HEADER_INCLUDED
#define LAPACK_HEADER_INCLUDED

#ifndef LAPACK_GLOBAL
#define LAPACK_GLOBAL(lcname,UCNAME)  lcname##_
#endif

#endif


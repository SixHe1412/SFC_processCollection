#ifndef _WIN32_EX_H_
#define _WIN32_EX_H_

#ifdef _WIN32
typedef char int8_t;
typedef unsigned char uint8_t;
typedef short int16_t;
typedef unsigned short uint16_t;
typedef int int32_t;
typedef unsigned int uint32_t;
typedef long long int64_t;
typedef unsigned long long uint64_t;

typedef unsigned int ssize_t;
typedef int socklen_t;

#include <WinSock2.h>
#include <time.h>
#include <ws2tcpip.h>
#define close closesocket

#else
#include <stdint.h>
#include <sys/time.h>
#endif


#ifdef _WIN32

#ifndef va_copy
#define va_copy(dst, src) ((void)((dst) = (src)))
#endif

#if defined( _WIN32 ) && !defined( __cplusplus )
#define inline __inline
#endif

#define EINPROGRESS WSAEINPROGRESS
#define EHOSTUNREACH WSAEHOSTUNREACH
#define EADDRNOTAVAIL WSAEADDRNOTAVAIL
//#define snprintf sprintf_s
#define strcasecmp strcmp
#define strncasecmp _strnicmp
#define strerror_r(errorno, buf, len) strerror_s(buf, len, errorno)

#endif



#endif

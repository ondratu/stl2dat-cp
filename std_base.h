/*======================================================================*
 *   TITLE: Standard types and declaration
 *   -----------------------------------------------------------------
 *   Advice for porting:
 *   -------------------------
 *   Take great care. Not portable at all
 =======================================================================*/
#ifndef __STD_BASE_H
#define __STD_BASE_H

#include <string.h>
#include <stdio.h>

/*-
 ! Following two define for declaring exprted functions
 ! in a way compatible with the job 'int'
 */
#define int_auto                  
#define DEF_VAL(a)                  

/*-
 ! Different system include
 */
#define DECL_SELF_CONSTRUCTOR(class_name)
#define SELF_CONSTRUCTOR(class_name)


#ifdef WIN32
#define STATIC_DOS
#define DOSFILE
// __STDC__ indicate ANSI compatibility,
// MSTDC also, but is used only by MTEL source code
#define MSTDC
#ifdef WINNT
#define MFAR
#define USEDBYDLL __declspec(dllimport)
#define DLL_EXPORT_VAR __declspec(dllexport)
#define DLL_EXPORT __declspec(dllexport)
// Used by Acis
#define NT
#define MSVC
#define _X86_
//#define _WINDOWS
#else
#define BIT16
#define MFAR __far
#define USEDBYDLL
#define DLL_EXPORT_VAR
#define DLL_EXPORT
#endif	// WINNT

#include <stdlib.h>
#include <fstream> // Define ifstream and ofstream
using namespace std;
#include <iomanip>
extern char* form( char *format, ... );
inline long sqr(long x) {return(x*x);}
#else
#define STATIC_DOS static
#define MFAR
#define USEDBYDLL
#define DLL_EXPORT_VAR
#define DLL_EXPORT
#endif

/*-
 ! Declaring standard types and values for all computers
 ! see chap L.2
 */
#define i32_to_i16(i) i

typedef bool t_Mbool;
#ifdef WIN32
#pragma warning(disable: 4237)
#pragma warning(disable: 4273)
#pragma warning(disable: 4786)
#endif
#ifndef VC50
//const t_Mbool true = 1;
//const t_Mbool false = 0;
#endif

typedef unsigned char t_Mbyte;
typedef short i16;
typedef long i32;
typedef unsigned long ulong;

typedef int CStringuct_id;

/*-
 !    Real are double. It is defined here for str.h
 */
typedef double t_Mreal;           

typedef t_Mbyte byte2;

typedef void *univ_ptr;
typedef i32 (*univ_ptr_fct)(); // Pointeur de fonction sans parametres

/*-
 ! defined in set.h
 */
class set32;
class set256;
class set1024;

typedef int t_index;
static const univ_ptr NIL =0;
#define IS_NIL(PRM) ((PRM)==0)


/*-
 ! Transaction type. see vmm.h for use ( it is rare )
 */
/*-
 ! 0 value for database pointer. May be != 0 for some computers
 ! Use vmm1.h to manipulate t_Mp
 */

//#include <str.h>
/*-
 ! Caracter strings
 */

/*-
 ! About assert. Skip to assert end
 */
/* Define assert whitch terminate by a kill -11 instead of exit(1) */
#include <signal.h>

//extern CString tb_where;
#ifdef WIN32
#ifdef NDEBUG

#define assert(exp) ((void)0)

#else 
extern "C" {
_CRTIMP void __cdecl _assert(void *, void *, unsigned);
}

#define assert(exp) \
    ( (exp) ? (void) 0 : _assert("", __FILE__, __LINE__) )

#endif 
#else
#include <assert.h>
#endif

extern DLL_EXPORT void errlog_traceback();
extern DLL_EXPORT void errlog_traceback(char*, char*, i32);

#ifdef RELEASE
#define precondition( test) 
#define postcondition( test) 
#define invariant( test) 
#else
#if defined SUN || defined MSDOS
#define precondition(test) { if ( ! (test)) { errlog_traceback("Precondition", __FILE__, __LINE__); assert( false); }}
#define postcondition(test) { if ( ! (test)) { errlog_traceback("Postcondition", __FILE__, __LINE__); assert( false); }}
#define invariant(test) { if ( ! (test)) { errlog_traceback("Invariant", __FILE__, __LINE__); assert( false); }}
#else
#define precondition( test) assert( test)
#define postcondition( test) assert( test)
#define invariant( test) assert( test)
#endif
#endif

/*-
 ! end assert.
 */


#endif



#ifndef __inp_H__
#define __inp_H__


#include <stdio.h>
#include <stdlib.h>

/* Nonzero if inp_handler callback should accept lineno parameter. */
#ifndef INP_HANDLER_LINENO
#define INP_HANDLER_LINENO 0
#endif

/* Typedef for prototype of handler function. */
#if INP_HANDLER_LINENO
typedef int (*inp_handler)(void* data, const char* section,
                           const char* name, const char* value,
                           int lineno);
#else
typedef int (*inp_handler)(void* data, const char* section,
                           const char* name, const char* value);
#endif

/* Typedef for prototype of fgets-style reader function. */
typedef char* (*inp_reader)(char* str, int num, void* stream);

/* Parse given inp-style file. May have [section]s, name=value pairs
   (whitespace stripped), and comments starting with ';' (semicolon). Section
   is "" if name=value pair parsed before any section heading. name:value
   pairs are also supported as a concession to Python's configparser.

   For each name=value pair parsed, call handler function with given data
   pointer as well as section, name, and value (data only valid for duration
   of handler call). Handler should return nonzero on success, zero on error.

   Returns 0 on success, line number of first error on parse error (doesn't
   stop on first error), -1 on file open error, or -2 on memory allocation
   error (only when inp_USE_STACK is zero).
*/
int inp_parse(const char* filename, inp_handler handler, void* data);

/* Same as inp_parse(), but takes a FILE* instead of filename. This doesn't
   close the file when it's finpshed -- the caller must do that. */
int inp_parse_file(FILE* file, inp_handler handler, void* data);

/* Same as inp_parse(), but takes an inp_reader function pointer instead of
   filename. Used for implementing custom or string-based I/O (see also
   inp_parse_string). */
int inp_parse_stream(inp_reader reader, void* stream, inp_handler handler,
                     void* data);

/* Same as inp_parse(), but takes a zero-terminated string with the inp data
instead of a file. Useful for parsing inp data from a network socket or
already in memory. */
int inp_parse_string(const char* string, inp_handler handler, void* data);



/* Chars that begin a start-of-line comment. Per Python configparser, allow
   both ; and # comments at the start of a line by default. */
#ifndef INP_START_COMMENT_PREFIXES
#define INP_START_COMMENT_PREFIXES "#"
#endif

/* Nonzero to allow inline comments (with valid inline comment characters
   specified by inp_INLINE_COMMENT_PREFIXES). Set to 0 to turn off and match
   Python 3.2+ configparser behaviour. */
#ifndef INP_ALLOW_INLINE_COMMENTS
#define INP_ALLOW_INLINE_COMMENTS 1
#endif
#ifndef INP_INLINE_COMMENT_PREFIXES
#define INP_INLINE_COMMENT_PREFIXES "#"
#endif


int split (const char *str, char c, char ***arr);


#define INP_MAX_LINE 200

#define INP_inpTIAL_ALLOC 200




#ifdef __cplusplus
}
#endif

#endif /* __inp_H__ */

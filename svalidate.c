/** @file 
 * \brief a small set of validation functions for strings 
 * \author Philip Papadopoulos
 * \email ppapadopoulos@ucsd.edu
 * \date 11 Jan 2018
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "svalidate.h"

/**
 * \fn char * trim (char * input)
 * \param input 
 * \brief trim leading and trailing white space from input string
 * \return pointer to beginning of trimmed string
 * \note Potential side-effect modifies:  can modify the input to place c-terminator 
 *       in the appropriate location to trim trailing whitespace
 */
char * trim (char * input)
{
	/* check null string */
	if (!input) return input;

	char * rstr = input;
	char * p = input;
	/* trim leading space */
	while ( isspace(*p)) p++;
	rstr = p;
	
	/* trim trailing white space */
	p = rstr + strlen(rstr) - 1;
	while ( isspace(*p)) p--;
	*(p+1) = '\0';
	return rstr;
}
/**
 * \fn char * trimCopy( char *input)
 * \param input 
 * \brief trim leading and trailing white space from input string
 * \return a malloced copy of trimmed string
 */
char * trimCopy( char *input)
{
	/* check null string */
	if ( !input ) return input;
	char * rstr;
	while ( isspace(*input) ) input++;
	if ( strlen(input) > 0 )
	{
		rstr = malloc((strlen(input) + 1)*sizeof(char));
		strcpy(rstr,input);
		return trim(rstr);
	}
	return (char *) NULL;
}
	
/**
 * \fn int isInteger(char * input)
 * \param input - string that is trimmed of leading/trailing whitespace 
 * \brief determine if input could be converted to an int 
 * \return 1 if string could be converted to an int, 0 otherwise 
 */
int isInteger(char * input)
{
	/* check null string */
	if (!input) return FALSE;
	if (*input == '+' || *input == '-') input++;
	do {
		if ( ! isdigit(*input)) return FALSE;
	}
	while (*(++input));
	return TRUE;
	
}
/**
 * \fn int isFloat(char * input)
 * \param input - string that is trimmed of leading/trailing whitespace 
 * \brief determine if input could be converted to a float
 * \return 1 if string could be converted to a float, 0 otherwise 
 */
int isFloat(char * input)
{
	float tval;
	char *endp;
	tval = strtof(input,&endp);
	// Check for various errors (see man strtof) 
	if ( *endp != '\0' ) return FALSE;
	if ( tval == 0 && endp == input ) return FALSE;
	if ( errno == ERANGE ) return FALSE;

	return TRUE;
}
/**
 *
 * \fn int isDouble(char * input)
 * \param input - string that is trimmed of leading/trailing whitespace
 * \brief determine if input could be converted to a double
 * \return 1 if string could be converted to a double, 0 otherwise
 */
int isDouble(char * input)
{
	float tval;
	char *endp;
	tval = strtod(input,&endp);
	// Check for various errors (see ma strtod)
	if ( *endp != '\0' ) return FALSE;
	if ( tval == 0 && endp == input ) return FALSE;
	if ( errno == ERANGE) return FALSE;

	return TRUE;
}

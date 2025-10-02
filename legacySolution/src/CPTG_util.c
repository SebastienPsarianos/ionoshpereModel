/* Use the preprocessor to include CPTG header file CPTG_util.h. */

#include "CPTG_util.h"

/*******************************************************************

 CPTG_read_line: This subroutine reads in a line of text from the
                 standard input file or terminal.

 Variable description:

 (call statement parameter list)

 buffer              Line of text.

 Begin subroutine CPTG_read_line. */

void CPTG_read_line(buffer)

/* Variable declaration. */

char buffer[];
{

  /* Local variable declaration. */

  int i;
  char single_char;

  /* Get the character string using getchar. */

  i=0;
  do
    {
       single_char = getchar();
       buffer[i] = single_char;
       ++i;
    }
  while(single_char != '\n');

  buffer[i-1] = '\0';

  return;

/* End subroutine CPTG_read_line. */

}

/*******************************************************************

 CPTG_concatenate: This subroutine concatenates two character strings.

 Variable description:

 (call statement parameter list)

 string1, string2          Input character strings.

 result                    Result of concatenation.

 Begin subroutine CPTG_concatenate. */

void CPTG_concatenate(string1, string2, result)

/* Variable declaration. */

char string1[], string2[], result[];
{

  /* Local variable declaration. */

  int i, j;

  /* Copy string1 to result. */

  for ( i = 0; string1[i] != '\0'; ++i )
    result[i] = string1[i];

  /* Copy string2 to result. */

  for ( j = 0; string2[j] != '\0'; ++j )
    result[i+j] = string2[j];

  /* Terminate the concatenated string with a null. */

  result[i+j] = '\0';

  return;

/* End subroutine CPTG_concatenate. */

}

/*******************************************************************

 CPTG_copy: This subroutine copies one character string to another.

 Variable description:

 (call statement parameter list)

 string_from, string_to          Character strings.

 Begin subroutine CPTG_copy. */

void CPTG_copy(string_from, string_to)

/* Variable declaration. */

char string_from[], string_to[];
{

  /* Local variable declaration. */

  int i;

  /* Copy string_from to string_to. */

  for ( i = 0; string_from[i] != '\0'; ++i )
    string_to[i] = string_from[i];

  string_to[i] = '\0';

  return;

/* End subroutine CPTG_copy. */

}

/*******************************************************************

 CPTG_compare: This function compares two character strings
               and determines whether they are equal or not,
               returning 1 (TRUE) if in fact the two strings
               are identical and 0 (FALSE) if they are not.

 Variable description:

 (call statement parameter list)

 string1, string2          Input character strings.

 Begin subroutine CPTG_compare. */

int CPTG_compare(string1, string2)

/* Variable declaration. */

char string1[], string2[];
{

  /* Local variable declaration. */

  int i=0, answer;

  /* Compare string1 and string2 */

  while ( string1[i] == string2[i] &&
          string1[i] != '\0' && string2[i] != '\0' )
     ++i;

  if (string1[i] == '\0' && string2[i] == '\0') {
     answer = 1;
  } else {
     answer = 0;
  } /* endif */

  return(answer);

/* End function CPTG_compare. */

}

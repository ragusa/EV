/**
 * \file Utilities.h
 * \brief Provides miscellaneous utilities.
 */
#ifndef Utilities_h
#define Utilities_h

#include <iostream>
#include <string>

namespace utilities
{
/**
 * \brief Wraps a string after a maximum number of characters per line.
 *
 * \param[in] s  string
 */
/*
std::string wrap_line(std::string s)
{
  // maximum number of characters per line
  const unsigned int line_width = 80;

  for (unsigned int i = 1; i <= s.length(); ++i)
  {
    // get character
    char c = s[i-1];

    int space_count = 0;

    // if line width has been met
    if ((i % line_width) == 0)
    {
      // if the character is not a space
      if (c != ' ')
      {
        // backtrack until space before word is found
        for (int j = i-1; j > -1; j--)
        {
          // if space before word is found
          if (s[j] == ' ')
          {
            // insert a newline
            s.insert(j, space_count, '\n');
            break;
          }
          else space_count++;
        }
      }
    }
  }

  return s;
}
*/

/**
 * \brief Issues a warning in a warning color.
 *
 * \param[in] warning  warning message
 */
void issue_warning(const std::string & warning)
{
  std::cout << "\x1b[33m"
            << "WARNING:" << std::endl
            << warning << "\x1b[0m" << std::endl;
}
}

#endif

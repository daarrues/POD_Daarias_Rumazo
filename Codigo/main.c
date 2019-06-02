//------------------------------------------------------------------------------
//                              Main
//------------------------------------------------------------------------------
// POD: Preliminary Orbit Determination
/**
 * @file main.c
 * @author Daniel Arias Ruiz-Esquide y Rubén Mazo Tomás
 * @date Created: 2019/06/02
 *
 * Este fichero contiene la función principal con el menú para ejecutar los
 * ejemplos de M. Mahooti y los test de integración.
 */
//------------------------------------------------------------------------------
#include "PreliminaryOrbitDetermination/examples.h"
#include "Test/test.h"
#include <stdio.h>

int main()
{
  test();
  printf("Example 5:\n");
  example5();
  return 0;
}

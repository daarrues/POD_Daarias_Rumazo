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

//------------------------------------------------------------------------------
//  int main()
//------------------------------------------------------------------------------
/**
 * Función principal, puede ejecutar los ejemplos o test.
 * Pregunta por la opción a ejecutar:
 * 1, 2, 3, 5, 6 y 7: ejemplo correspondiente.
 * 9: test (todos).
 * 0: salir.
 */
//------------------------------------------------------------------------------
int main()
{
  int op = 0;
  do
  {
    printf("---[MENU]---\n");
    printf("1. Example 1\n");
    printf("2. Example 2\n");
    printf("3. Example 3\n");
    printf("5. Example 5\n");
    printf("6. Example 6\n");
    printf("7. Example 7\n");
    printf("\n9.  Test\n");
    printf("\n0.  Exit\n");
    printf("------------\n");
    printf("Option: ");
    scanf("%d", &op);
    printf("------------\n");
    switch(op)
    {
      // Examples
      case 1:
        printf("Example 1:\n");
        example1();
        break;
      case 2:
        printf("Example 2:\n");
        example2();
        break;
      case 3:
        printf("Example 3:\n");
        example3();
        break;
      case 5:
        printf("Example 5:\n");
        example5();
        break;
      case 6:
        printf("Example 6:\n");
        example6();
        break;
      case 7:
        printf("Example 7:\n");
        example7();
        break;
      case 9:
        test();
        break;
      case 0:
        printf("Exiting...\n");
        return 0;
      default:
        printf("ERROR: %d is not a valid option!\n", op);
        break;
    }
    printf("\n");
  } while (1);
}

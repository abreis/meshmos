//

#include "object.h"
#include "simulator.h"

#include <iostream>


MOSstat::MOSstat () {

}

int MOSstat::command(int argc, const char*const* argv) {

	// do this
	if ( argc == 1 && strcmp (argv[0], "on") == 0 ) {
		return 0;  // TCL_OK

	// do that
	} else if ( argc == 1 && strcmp (argv[0], "off") == 0 ) {

		return 0;  // TCL_OK
	}


	fprintf (stderr, "invalid mosstat command: '");
	for ( int i = 0 ; i < argc ; i++ ) fprintf (stderr, "%s ", argv[i]);
	fprintf (stderr, "'\n");

	return 1; // TCL_ERROR
}

#include <stdio.h>
#include <string.h>

#include "interface_artic_c.h"


int main()
{
	printf("running example....\nartic_c help = ");

	enum { kMaxArgs = 64 };
	int argc = 0;
	char *argv[kMaxArgs];

	 char commandLine[250] = "artic_c -h";
	char *p2 = strtok(commandLine, " ");
	
	while (p2 && argc < kMaxArgs-1)
	  {
	    argv[argc++] = p2;
	    p2 = strtok(0, " ");
	  }
	argv[argc] = 0;
	
	init_artic(argc,argv);

	return 0;
}


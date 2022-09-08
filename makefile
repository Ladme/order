all: order.c ordermap.c general.c
	make order groan=$(groan)
	make ordermap groan=$(groan)

order: order.c general.c
	gcc order.c general.c -I$(groan) -L$(groan) -D_POSIX_C_SOURCE=200809L -o order -lgroan -lm -std=c99 -pedantic -Wall -Wextra -O3 -march=native

ordermap: ordermap.c general.c
	gcc ordermap.c general.c -I$(groan) -L$(groan) -D_POSIX_C_SOURCE=200809L -o ordermap -lgroan -lm -std=c99 -pedantic -Wall -Wextra -O3 -march=native

install:
	if [ -f order ]; then cp order ${HOME}/.local/bin; fi
	if [ -f ordermap ]; then cp ordermap ${HOME}/.local/bin; fi

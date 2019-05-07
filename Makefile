CFLAGS += -O3 -Wno-implicit-int -Wno-return-type \
	-Wno-implicit-function-declaration -Wno-format
LIBS = -lm -lX11
DESTBIN = /usr/local/bin
DESTMAN = /usr/local/share/man/man1

xmartin: xmartin.c Makefile
	$(CC) $(CFLAGS) -o xmartin xmartin.c $(LIBS)

install: xmartin
	install -m 755 xmartin $(DESTBIN)
	install -m 755 xmartin+ $(DESTBIN)
	install xmartin.1 $(DESTMAN)
      
clean:
	$(RM) xmartin

uninstall:
	$(RM) $(DESTBIN)/xmartin $(DESTBIN)/xmartin+ $(DESTMAN)/xmartin.1

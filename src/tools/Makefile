
all:
	( $(MAKE) build_registry CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" )
	( $(MAKE) build_input_gen CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" )
	( $(MAKE) build_string_consts CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" )

build_input_gen:
ifdef MPAS_TOOL_DIR
	@echo "*** Using MPAS tools from ${MPAS_TOOL_DIR} ***"
	(cp $(MPAS_TOOL_DIR)/input_gen/namelist_gen input_gen/)
	(cp $(MPAS_TOOL_DIR)/input_gen/streams_gen input_gen/)
else
	@echo "*** Building MPAS tools from source ***"
	(cd input_gen; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)")
endif

build_registry:
ifdef MPAS_TOOL_DIR
	@echo "*** Using MPAS tools from ${MPAS_TOOL_DIR} ***"
	(cp $(MPAS_TOOL_DIR)/registry/parse registry/)
else
	@echo "*** Building MPAS tools from source ***"
	(cd registry; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)")
endif

build_string_consts:
ifdef MPAS_TOOL_DIR
	@echo "*** Using MPAS tools from ${MPAS_TOOL_DIR} ***"
	(cp $(MPAS_TOOL_DIR)/string_consts/gen_string_consts string_consts/)
else
	@echo "*** Building MPAS tools from source ***"
	(cd string_consts; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)")
endif

clean:
	(cd input_gen; $(MAKE) clean)
	(cd registry; $(MAKE) clean)
	(cd string_consts; $(MAKE) clean)

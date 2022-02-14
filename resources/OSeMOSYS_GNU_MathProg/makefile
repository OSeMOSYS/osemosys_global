VERSION := $(shell git describe)
BUILDDIR := osemosys_gnu_mathprog_$(VERSION)

release: $(BUILDDIR).zip

$(BUILDDIR).zip: model_files data_files scripts $(BUILDDIR)/README.html $(BUILDDIR)/LICENSE
	zip -r $(BUILDDIR).zip $(BUILDDIR)

model_files: $(BUILDDIR)/osemosys.txt $(BUILDDIR)/osemosys_short.txt $(BUILDDIR)/osemosys_fast.txt

data_files:	$(BUILDDIR)/simplicity.txt

scripts: $(BUILDDIR)/scripts

$(BUILDDIR)/scripts: $(BUILDDIR)
	cp -R scripts $(BUILDDIR)/scripts

$(BUILDDIR)/simplicity.txt:	$(BUILDDIR)
	cp tests/simplicity.txt $(BUILDDIR)/simplicity.txt

$(BUILDDIR)/osemosys.txt: $(BUILDDIR)
	cp src/osemosys.txt $(BUILDDIR)/osemosys.txt

$(BUILDDIR)/osemosys_short.txt: $(BUILDDIR)
	cp src/osemosys_short.txt $(BUILDDIR)/osemosys_short.txt

$(BUILDDIR)/osemosys_fast.txt: $(BUILDDIR)
	cp src/osemosys_fast.txt $(BUILDDIR)/osemosys_fast.txt

$(BUILDDIR)/README.html: $(BUILDDIR)
	pandoc -f markdown -t html src/README.md -o $(BUILDDIR)/README.html

$(BUILDDIR)/LICENSE: $(BUILDDIR)
	cp LICENSE $(BUILDDIR)/LICENSE

$(BUILDDIR):
	mkdir -p "$(BUILDDIR)"

.PHONY: clean release all
all: release
clean:
	rm -f $(BUILDDIR).zip
	rm -rf $(BUILDDIR)
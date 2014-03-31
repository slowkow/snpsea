HTMLOPT = --standalone --toc --toc-depth=2 --template skeleton --mathjax

.PHONY: all manual clean

all: manual


manual: SNPsea_manual.md
	perl -pe 's/\.pdf/.png/g' $? | pandoc $(HTMLOPT) -o SNPsea_manual.html


clean:
	rm -f SNPsea_manual.html

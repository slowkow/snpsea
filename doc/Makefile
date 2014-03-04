HTMLOPT = --standalone --toc --toc-depth=2 --template skeleton --mathjax
PDFOPT = --standalone --normalize --toc --toc-depth=2
VARS += --variable geometry:margin=1in
VARS += --variable fontsize:11pt

.PHONY: all manual clean

all: manual


manual: SNPsea_manual.md
	perl -pe 's/\.pdf/.png/g' $? | pandoc $(HTMLOPT) -o SNPsea_manual.html
	pandoc $(PDFOPT) $(VARS) -i $? -o SNPsea_manual.pdf


clean:
	rm -f SNPsea_manual.pdf SNPsea_manual.html

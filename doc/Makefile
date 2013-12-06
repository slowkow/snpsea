HTMLOPT = -s --toc --toc-depth=2 --template nice
PDFOPT = -s --normalize --toc --toc-depth=2
VARS += --variable geometry:margin=1in
VARS += --variable fontsize:11pt

.PHONY: all supplement manual clean

all: supplement manual


supplement: SNPsea_supplement.md
	pandoc $(PDFOPT) $(VARS) $? -o SNPsea_supplement.pdf


manual: SNPsea_manual.md
	perl -pe 's/\.pdf/.png/g' $? | pandoc $(HTMLOPT) -o SNPsea_manual.html
	pandoc $(PDFOPT) $(VARS) -i $? -o SNPsea_manual.pdf


clean:
	rm -f SNPsea_manual.pdf SNPsea_manual.html
	rm -f SNPsea_supplement.pdf
